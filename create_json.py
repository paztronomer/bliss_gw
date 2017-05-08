"""Simple script that creates JSON files for CTIO 4m, given RA,DEC and
whether is half night or entire one.
"""
import os
import sys
import time
import subprocess
import shlex
import argparse
import numpy as np
import astropy
import astropy.units as apy_u
import astropy.coordinates as apy_coord
import astropy.time as apy_time
import scipy.signal as signal #v0.19
import scipy.interpolate as interpolate
#others?
#import json
import logging

"""STEPS
1) CTIO coordinates --DONE
2) set the night window for some half or full night --DONE
3) restrict the window with the available dome pointing options --DONE
4) open text file with ra,dec,nite,first/second/full
5) decide the order of the exposures
6) populate JSON dictionary
7) write out the JSON file
"""

class Toolbox():
    @classmethod
    def lsq_interp(cls,x,y,degree=4):
        """Caller fot the scipy least squares method interpolator
        To call scipy.interpolate_make_lsq_spline():
        - x: abcissas
        - y: ordinates
        - t: knots, array-like with shape(n+k+1)
        - k: b-spline degree
        - w: weights for spline fitting
        - axis: interpolation axis, default is 0
        - check_finite: whether to check if the inputs contains only finite
        elements
        NOTES:
        (*) number of data points must be larger than the spline degree
        (*) knots must be a subset of data points of x[j] such that
        t[j] < x[j] < t[j+k+1], for j=0,1..n-k-2
        (*) degree 4 works slightly better than lower values
        """
        naux = np.int(np.ceil(x.shape[0]*0.1))
        #by hand I set the use of all but 2 points at the ends
        p1,p2 = 2,x.shape[0]-3
        t = x[p1:p2:naux]
        #using np.r_[] setup an array
        t = np.r_[(x[0],)*(degree+1),t,(x[-1],)*(degree+1)]
        lsq_spline = interpolate.make_lsq_spline(x,y,t,degree)
        return lsq_spline

    @classmethod
    def lsq_spline1D(cls,x,y):
        from scipy.interpolate import LSQUnivariateSpline, UnivariateSpline
        naux = np.int(np.ceil(x.shape[0]*0.1))
        p1,p2 = 5,x.shape[0]-4
        #t = x[p1:p2:naux]
        t = x[::10]
        degree = 3
        t = np.r_[(x[0],)*(degree+1),t,(x[-1],)*(degree+1)]
        spl = LSQUnivariateSpline(x,y,t)
        return spl

    @classmethod
    def delta_hr(cls,time_arr):
        """Method to return a round amount of hours for a astropy TimeDelta
        object, calculated from the peak-to-peak value given an astropy
        Time array
        Inputs:
        - time_arr: arreay of Time objects
        Returns
        - N: integer representing the round value of the number of hours the
        time interval contains.
        - M: integer representing the rounf dumber of minutes the interval contains
        """
        N = np.round(time_arr.ptp().sec/3600.).astype(int)
        M = np.round(time_arr.ptp().sec/60.).astype(int)
        if (np.abs(N) < 1) and (np.abs(N) >= 0):
            logging.warning("Observing window is smaller than 1 hour")
        elif N > 24:
            logging.warning("Observing window is longer than 1 day")
        elif N < 0:
            logging.error("Error: time range is negative")
            exit(1)
        return (N,M)


class Telescope():
    @classmethod
    def site(cls,name=None):
        #method from: http://docs.astropy.org/en/stable/api/astropy.coordinates.
        #EarthLocation.html
        #geodetic units
        #- lat/long: [deg,min,sec]
        #- altitude: mt
        #- ellipsoid: WGS84 where South and West are negative

        #CTIO geodetic
        #data from: http://www.ctio.noao.edu/noao/node/2085
        geod_lat =  [-30,10,10.78]
        geod_long = [-70,48,23.49]
        geod_h = 2241.4
        geod_ellip = "WGS84"

        if name is not None:
            coo = apy_coord.EarthLocation.of_site(name)
        else:
            sign_lat,sign_long = np.sign(geod_lat[0]),np.sign(geod_long[0])
            geod_lat = [abs(x)/np.float(60**i) for i,x in enumerate(geod_lat)]
            geod_long = [abs(x)/np.float(60**i) for i,x in enumerate(geod_long)]
            geod_lat = sign_lat*sum(geod_lat)
            geod_long = sign_long*sum(geod_long)
            coo = apy_coord.EarthLocation.from_geodetic(
                geod_long,geod_lat,height=geod_h,ellipsoid=geod_ellip)
        return coo

    @classmethod
    def horizon_limits_tcs(cls,ra_zenith,time_interp):
        """Method to return an interpolator, based on the coordinates given
        by CTIO as the horizon limits. The coordinates from CTIO are in
        HOUR_ANGLE,DEC so the HOUR_ANGLE must be added/substracted to the RA
        of the zenith at a particular time, to get the real limits.
        HOUR_ANGLE = 0 == RA at ZENITH
        Returns:
        - an interpolator object, class scipy.interpolate._bsplines.BSpline
        When this interpolator is applied, returns an array
        """
        HAngle= [5.25,5.25,5.25,5.25,5.25,5.25,5.25,5.25,5.25,5.25,5.25,
                5.25,5.25,5.12,4.96,4.79,4.61,4.42,4.21,3.98,3.72,3.43,3.08,
                2.64,2.06,1.10,0.00]
        Dec_d = [-89.00,-85.00,-80.00,-75.00,-70.00,-65.00,-60.00,-55.00,
                -50.00,-45.00,-40.00,-35.00,-30.00,-25.00,-20.00,-15.00,
                -10.00,-05.00,00.00,05.00,10.00,15.00,20.00,25.00,30.00,
                35.00,37.00]
        alt_d = [30.4,31.0,31.6,31.9,32.0,31.8,31.3,30.6,29.6,28.3,26.9,
                25.2,23.4,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,
                23.0,23.0,23.0,23.0,23.0]
        HAngle,Dec_d = np.array(HAngle),np.array(Dec_d)
        HAngle = apy_coord.Angle(HAngle,unit=apy_u.h)
        #interpolate initial grid
        #lsq_spl = Toolbox.lsq_interp(np.array(dec_d),np.array(hra_h))

        #Steps:
        #1) for each of the time stamps, transform HourAngle to RA!!!
        #2) return the interpolator, using RA instead of
        #Note: I will calculate only the upper envelope, as this is symmetric
        shift = []
        for zen in ra_zenith[:,0]:
            tmp_ra = HAngle.degree + zen
            lsq_dec = Toolbox.lsq_interp(np.array(Dec_d),np.array(tmp_ra))
            shift.append((lsq_dec,tmp_ra,zen))
        #Usage: xs, tmp_lsq(xs)
        return shift

    @classmethod
    def horizon_limits_plc(cls):
        #limits: [ha_west,ha_east,dec_south,dec_north]
        ha_d = [84.25,-89.25,0.00,0.00]
        dec_d = [-30.0,-22.0,-95.14,40.31]
        alt_d = [245.00,109.00,180.00,0.00]
        az_d = [18.97,11.40,24.86,19.68]
        return False

    @classmethod
    def zenith_ra(cls,time_arr,site):
        """Calculates the RA of the azimuth, at a given time, at a given site.
        Time must be given in UTC.
        Inputs:
        - array of times (astropy.time class) containing the set of time stamps
        between the borders. The shape is (interpolation_inside_interval,1)
        - site: location of the observing site
        Returns:
        - array of same dimensions as the input array. Each element is a
        float for the RA of the zenith at given time and site, in degrees.
        """
        #create the azimuth coordinates and then use them to get the RA
        g = lambda x: apy_coord.SkyCoord(alt=90.*apy_u.deg,az=0.*apy_u.deg,
                                        obstime=x,location=site,
                                        frame="altaz").icrs.ra.degree
        h = lambda x: np.array(map(g,x))
        res = np.array(map(h,time_arr))
        return res


class Schedule():
    @classmethod
    def eff_night(cls,day_ini,earth_loc):
        """This method calculates the time range between the Sun being
        at -14deg below the horizon, for a single night of observation
        Inputs:
        - day_ini: noon of the initial day of observation
        - earth_loc: location of the observing site
        Returns:
        - array of astropy.time.core.Time elements, containing begin,
        middle,and end of the night
        """
        #Sun position is in GCRS
        aux = day_ini + np.linspace(0,24,1440) * apy_u.hour
        sun1 = apy_coord.get_sun(aux)
        altaz_frame = apy_coord.AltAz(obstime=aux,location=earth_loc)
        sun1_altaz = sun1.transform_to(altaz_frame)
        #transform AltAz to degrees. Use vectorization
        vs,v1,v2,v3 = apy_coord.Angle(sun1_altaz.alt).signed_dms
        todeg = lambda w,x,y,z: w*(x + y/np.float(60) + z/np.float(60**2))
        vect = np.vectorize(todeg,otypes=["f4"])
        sun1_deg = vect(vs,v1,v2,v3)
        #Sun position closest to -14deg, use local minima
        idx_minim = signal.argrelextrema(np.abs(sun1_deg + 14),np.less)
        if idx_minim[0].shape[0] != 2:
            logger.error("Error setting Sun at -14deg")
            exit(1)
        #save time for begin, middle, and end of the night
        ntime = aux[idx_minim]
        midnight = apy_time.Time(apy_time.Time(np.mean(ntime.jd),format="jd"),
                                scale="utc",format="iso")
        ntime = np.insert(ntime,1,midnight)
        return ntime

    @classmethod
    def scan_night(cls,time_kn,Nstep=24):
        """Use the 3 times for begin,middle,end of the night to create a
        set of intermediate values.
        Inputs:
        - time_kn = 1D array with 2 or 3 astropy Time entries, representing
        begin, (optional: middle), end of the night
        Returns:
        - array of shape (M,N) N=1 being number of section of the night
        and M the number of time intervals the section want to be divided for a
        later interpolation. Each time is an object of the class
        astropy.time.core.Time
        """
        #to iso format
        ft = lambda x: apy_time.Time(apy_time.Time(x,format="jd"),
                                    scale="utc",format="iso")#.iso
        #cases for 1 and 2 divisions of the night
        tjd = map(lambda x: x.jd, time_kn)
        if (time_kn.shape[0] in [2,3]):
            tjd_uni = np.linspace(np.min(tjd),np.max(tjd),Nstep+1)
            iso_uni = np.array(map(ft,tjd_uni))#,dtype=np.dtype("|S50"))
            res = iso_uni
        else:
            logging.error("Scan method needs 2-3 astropy-Time inputs")
            exit(1)
        res = res[:,np.newaxis]
        return res

    @classmethod
    def avoid_moon(cls,day_ini,earth_loc):
        #set a radius of N degrees from moon, if full??
        moon1 = apy_coord.get_moon(day_ini,earth_loc)
        return False

    @classmethod
    def point(cls,site_name=None,
            utc_minus_local=3,
            begin_day="2017-04-13 12:00:00",
            obs_interval="full"):
        """This method has few main steps:
        1) gets the site location from the Telescope class
        2) gets the range of hours every night lasts (to get the splitting for
        half observation nights)
        3) pick the coordinates of one field and gets its AltAz position at
        a given time
        PENDING:
        - if the field is visible, get the airmass (use threshold by user)
        and select the order of observing
        Note:
        - time from user must be feed in local time, but processing is made in
        UTC. The results are given in UTC time by astropy
        """
        #Earth location
        site = Telescope.site(name=site_name)
        #some SkyCoordinate, from Event 3
        ra,dec = [0.430246,-35.4841]
        #the RA DEC coordinates by convenience must be translated to
        #altitude(0-90) azimuth(0-360,N->E)

        #time range: must be set using the Sun elevation
        #now I"ll use april 13, from 8.34pm to 6.57am
        #Chile: UTC-3; Daylight save: UTC-4 -> local = UTC - <3,4>
        #deltaUTC = utc_diff*apy_u.hour
        #local_t1 = apy_time.Time("2017-04-13 20:34:00")
        #local_t2 = apy_time.Time("2017-04-14 06:57:00")
        #trange = [local_t1,local_t2]
        #trange = [x - deltaUTC for x in trange]

        """use UTC time plus the local site to get the obs window"""
        #starting at noon, scan for 24hrs searching for the Sun at -14deg,
        #returns an array on which time entries are astropy.time.core.Time
        deltaUTC = utc_minus_local*apy_u.hour
        t_utc = apy_time.Time(begin_day) + deltaUTC
        t_window = Schedule.eff_night(t_utc,site)
        if obs_interval == "first":
            t_window = t_window[:-1]
        elif obs_interval == "second":
            t_window = t_window[1:]
        elif (obs_interval == "full") or (obs_interval is None):
            pass
        else:
            logging.error("Interval must be: first, second, full")

        """inside the obs window, define a grid of times to be used as steps
        for the calculations of obs window"""
        #given the borders or the observation window for a specific night,
        #interpolate as many points as the number oef hours in the interval,
        #as defined by the input borders.
        #Returns an array of astropy.time.core.Time entries
        N_hr,N_min = Toolbox.delta_hr(t_window)
        t_interp = Schedule.scan_night(t_window,Nstep=N_hr)

        """calculate the RA of the zenith in each of the time stamps"""
        #for each of the time stamps, find the zenith RA. Returns an array of
        #same shape as the array of interpolated times
        zen_ra = Telescope.zenith_ra(t_interp,site)

        """giving RA for the zenith at diff times, locate the CTIO horizon"""
        #translate the borders found in CTIO PDF document to RA-DEC, initially
        #given in HourAngle,Dec
        fx_dec = Telescope.horizon_limits_tcs(zen_ra,t_interp)
        exit()


        #get the telescope horizon limits, in AltAz coordinates
        #lim_alt,lim_az,lim_altfix = Telescope.horizon_limits_tcs(t_window,site)

        #give time range with more entries? yes, because an object can
        #be at the observability window maybe for only a part of the night

        #transformation 1: non considering obstime in RADEC initialization
        radec = apy_coord.SkyCoord(ra=ra,dec=dec,frame="icrs",
                                unit=(apy_u.deg,apy_u.deg))
        f = lambda x : radec.transform_to(apy_coord.AltAz(obstime=x,
                                                        location=site))
        #list of objects AltAz
        altaz = [f(a) for a in trange]
        #to transform AltAz to angles in degrees, use Angle(altaz[0].alt).dms


        print apy_coord.Angle(altaz[0].alt).dms[:]
        print apy_coord.Angle(altaz[0].alt).is_within_bounds("0d", "360d")
        print altaz[0].alt
        print altaz[0].az
        print altaz[0].secz
        print apy_coord.Angle([-20, 150, 350]*apy_u.deg).is_within_bounds(
            -10*apy_u.deg,None)
        #transformation 2: considering obstime in RADEC
        #Nothing changes, so I will not use the additional arg





if __name__ == "__main__":

    print "starting..."
    Schedule.point()
    #Telescope.site()
