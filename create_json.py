'''Simple script that creates JSON files for CTIO 4m, given RA,DEC and
whether is half night or entire one.
'''
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

'''STEPS
1) CTIO coordinates --DONE
2) set the night window for some half or full night --DONE
3) restrict the window with the available dome pointing options
4) open text file with ra,dec,nite,first/second/full
5) decide the order of the exposures
6) populate JSON dictionary
7) write out the JSON file 
'''

class Toolbox():
    @classmethod
    def lsq_interp(cls,x,y,degree=3):
        '''Caller fot the scipy least squares method interpolator
        To call scipy.interpolate_make_lsq_spline():
        - x: abcissas
        - y: ordinates
        - t:knots, array-like with shape(n+k+1)
        - k: b-spline degree
        - w: weights for spline fitting
        - axis: interpolation axis, default is 0
        - check_finite: whetehr to check if the inputs contains only finite 
        elements
        NOTE: 
        (*) number of data points must be larger than the spline degree
        (*) knots must be a subset of data points of x[j] such that
        t[j] < x[j] < t[j+k+1], for j=0,1..n-k-2
        '''
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
        geod_ellip = 'WGS84'
        
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
    def horizon_limits_tcs(cls,time_range,earth_loc):
        '''The horizon limits must be calculated in a nightly basis, so this
        method receives the effective time range for the observing run. 
        The allowed coordinates are transformed to AltAz and returned
        '''
        #http://www.ctio.noao.edu/noao/content/Horizon-Limits
        #https://www.ctio.noao.edu/DocDB/0007/000717/002/Horizon%20Limits.pdf
        ra_hr = [5.25,5.25,5.25,5.25,5.25,5.25,5.25,5.25,5.25,5.25,5.25,
                5.25,5.25,5.12,4.96,4.79,4.61,4.42,4.21,3.98,3.72,3.43,3.08,
                2.64,2.06,1.10,0.00]
        dec_d = [-89.00,-85.00,-80.00,-75.00,-70.00,-65.00,-60.00,-55.00,
                -50.00,-45.00,-40.00,-35.00,-30.00,-25.00,-20.00,-15.00,
                -10.00,-05.00,00.00,05.00,10.00,15.00,20.00,25.00,30.00,
                35.00,37.00]
        alt_d = [30.4,31.0,31.6,31.9,32.0,31.8,31.3,30.6,29.6,28.3,26.9,
                25.2,23.4,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,23.0,
                23.0,23.0,23.0,23.0,23.0]
        min_alt = 23
        """
        ra_hr += map(lambda r: -1*r,ra_hr[:-1])
        dec_d += dec_d[:-1]
        alt_d += alt_d[:-1]
        #sorting to perform LSQ in dec,ra (inverted for simplicity)
        idx1 = np.argsort(dec_d)
        ra_hr = np.array(ra_hr)[idx1]
        dec_d = np.array(dec_d)[idx1]
        alt_d = np.array(alt_d)[idx1]
        """
        #interpolate
        lsq1 = Toolbox.lsq_interp(np.array(dec_d),np.array(ra_hr))
        dec_spl = np.linspace(min(dec_d),max(dec_d),1000)
        ra_spl = lsq1(dec_spl)
        #reflect data for negative ra
        dec_spl = np.concatenate([dec_spl,dec_spl[:-1]])
        ra_spl = np.concatenate([ra_spl,-1*ra_spl[:-1]])
        #sort in ra
        idx1 = np.argsort(ra_spl)
        ra_spl = np.array(ra_spl)[idx1]
        dec_spl = np.array(dec_spl)[idx1]
        #ra must be ]-5.25, 5.25[ deg
        #dec must be ]39,-89[ deg
        #elevation >23 deg
        #must use the combination of these 3 to get the horizon limit
        #translate them to the same frame: AltAz 
        '''alt_d = apy_coord.Angle(alt_d,unit=apy_u.deg)'''
        aux_radec = apy_coord.SkyCoord(ra=ra_spl,dec=dec_spl,frame='icrs',
                                unit=(apy_u.h,apy_u.deg))
        f = lambda x : aux_radec.transform_to(apy_coord.AltAz(obstime=x,
                                                        location=earth_loc))
        aux_altaz = [f(a) for a in time_range]
        #evaluate limits in each time border, returns arr [deg],[min],[sec]
        alt1 = [apy_coord.Angle(x.alt).signed_dms[:] for x in aux_altaz]
        az1 = [apy_coord.Angle(x.az).signed_dms[:] for x in aux_altaz]
        '''alt2 = [alt_d.signed_dms[:] for x in alt_d]'''

        #to transform to degrees
        sumdeg = lambda *w: w[0]*(w[1]+w[2]/np.float(60)+w[3]/np.float(60**2))
        vect = np.vectorize(sumdeg,otypes=['f8'])
        #for every time of the input time range
        alt1 = map(lambda x: vect(*x), alt1)
        az1 = map(lambda x: vect(*x), az1)

        #need to replicate in the initial altitude set of values, the number
        #of entries of the time range
        '''alt_d = [np.array(alt_d) for i in xrange(time_range.shape[0])]'''
      
        return alt1,az1,min_alt

        """
        '''MUST PERFORM THE LSQ ON RA DEC, BECAUSE IS EASIER AND FIXED!!
        '''
        exit()
        #sort ascending, using ra_hr
        idx = np.argsort(az1[2])[::]
        tmp_az1 = np.array(az1[2])[idx]
        tmp_alt1 = np.array(alt1[2])[idx]
        #interpolate
        #plt.plot(test_x, lsq_p(test_x), 'g-', lw=3, label='LSQ spline')
        #plt.plot(tmp_az1,tmp_alt1,'bo')
        #plt.show()
        #plt.plot(np.linspace(az1[0][0],az1[0][-1],1000), 
        #        lsq(np.linspace(az1[0][0],az[0][-1],1000)), 'o', 
        #        az1, alt1, '-')
        #plt.show()

        exit()

        xp = az1[2]
        fp = alt1[2]
        print np.interp(-10.,xp,fp)

        #interpolate
        import matplotlib.pyplot as plt
        from scipy import interpolate
        x = az1[2]
        y = alt1[2]
        f = interpolate.interp1d(x,y)
        xnew = np.arange(min(x),max(y),0.1)
        ynew = f(xnew)   # use interpolation function returned by `interp1d`
        plt.plot(x, y, 'o', xnew, ynew, '-')
        plt.show()
        exit()

        #aux_alt = [apy_coord.Angle(x.alt).signed_dms[:] for x in aux_altaz]
        #aux_az = [apy_coord.Angle(x.az).signed_dms[:] for x in aux_altaz]

        #transform AltAz to degrees. Use vectorization
        #a1s,a1v1,a1v2,a1v3 = apy_coord.Angle(sun1_altaz.alt).signed_dms
        #todeg = lambda w,x,y,z: w*(x + y/np.float(60) + z/np.float(60**2))
        #vect = np.vectorize(todeg,otypes=['f4'])
        #sun1_deg = vect(vs,v1,v2,v3)
        #print alt1
        #exit()
        
        import matplotlib.pyplot as plt
        for i in xrange(len(alt1)):
            plt.plot(az1[i],alt1[i],'bo')
            plt.plot(az1[i],alt_d[i],'ro')
        plt.show()
        
        exit()
        #to transform AltAz to angles in degrees, use Angle(altaz[0].alt).dms
        """
        """
        print apy_coord.Angle(altaz[0].alt).dms[:]
        print apy_coord.Angle(altaz[0].alt).is_within_bounds('0d', '360d')
        print altaz[0].alt
        print altaz[0].az
        print altaz[0].secz
        print apy_coord.Angle([-20, 150, 350]*apy_u.deg).is_within_bounds(
            -10*apy_u.deg,None)
        """

    @classmethod
    def horizon_blanco(cls):
        '''
        Overview:
        1) limits are given (CTIO webpage) in HourAngle,Dec. Dec is consistent
        with usual Declination. HourAngle is measured from RA at the zenith
        2) for each object with RA-DEC, use DEC as it is. For RA, must be 
        expressed in HourAngle, where HourAngle=0 is the Altitude of 90 deg
        3) having the object coordinates in RA-DEC expressed in HourAngle-DEC,
        must compare to the horizon limits given by CTIO
        4) after compare, cut by airmass
        '''
        pass

    @classmethod 
    def horizon_limits_plc(cls):
        #limits: [ha_west,ha_east,dec_south,dec_north]
        ha_d = [84.25,-89.25,0.00,0.00]
        dec_d = [-30.0,-22.0,-95.14,40.31]
        alt_d = [245.00,109.00,180.00,0.00]
        az_d = [18.97,11.40,24.86,19.68]
        return False

    @classmethod
    def other_restrictions():
        #https://www.ctio.noao.edu/DocDB/0007/000717/002/Horizon%20Limits.pdf
        return False
    
    @classmethod
    def zenith_ra(cls,time_arr,site):
        '''Calculates the RA of the azimuth, at a given time, at a given site.
        Time must be given in UTC.
        '''
        #create the azimuth coordinates and then use them to get the RA
        g = lambda x: apy_coord.SkyCoord(alt=90.*apy_u.deg,az=0.*apy_u.deg,
                                        obstime=x,location=site,
                                        frame='altaz').icrs.ra
        aux1 = map(g,time_arr)
        #print aux1[0]
        #print apy_coord.Angle(apy_coord.Angle(5.25,unit=apy_u.h),unit=apy_u.deg)
        #u1 = aux1[0] - apy_coord.Angle(5.25,unit=apy_u.h)
        return aux1


class Schedule():
    @classmethod 
    def eff_night(cls,day_ini,earth_loc):
        '''This method calculates the time range between the Sun being
        at -14deg below the horizon, for a single night of observation
        Inputs:
        - day_ini: noon of the initial day of observation
        - earth_loc: location of the observing site
        Returns:
        - array containing begin,middle,and end of the night
        '''
        #Sun position is in GCRS
        aux = day_ini + np.linspace(0,24,1440) * apy_u.hour
        sun1 = apy_coord.get_sun(aux)
        altaz_frame = apy_coord.AltAz(obstime=aux,location=earth_loc)
        sun1_altaz = sun1.transform_to(altaz_frame)
        #transform AltAz to degrees. Use vectorization
        vs,v1,v2,v3 = apy_coord.Angle(sun1_altaz.alt).signed_dms
        todeg = lambda w,x,y,z: w*(x + y/np.float(60) + z/np.float(60**2))
        vect = np.vectorize(todeg,otypes=['f4'])
        sun1_deg = vect(vs,v1,v2,v3)
        #Sun position closest to -14deg, use local minima
        idx_minim = signal.argrelextrema(np.abs(sun1_deg + 14),np.less)
        if idx_minim[0].shape[0] != 2: 
            logger.error('Error setting Sun at -14deg')
            exit(1)
        #save time for begin, middle, and end of the night
        ntime = aux[idx_minim]
        midnight = apy_time.Time(apy_time.Time(np.mean(ntime.jd),format='jd'),
                                scale='utc',format='iso')
        ntime = np.insert(ntime,1,midnight)
        return ntime
        
    @classmethod
    def scan_night(cls,time_kn,Nstep=500):
        '''Use the 3 times for begin,middle,end of the night to create a
        set of intermediate values
        Inputs:
        - time_kn = list of times objects, astropy Time
        '''
        #to iso format
        ft = lambda x: apy_time.Time(apy_time.Time(x,format='jd'),
                                    scale='utc',format='iso')
        if len(time_kn) == 2:
            tjd = map(lambda x: x.jd, time_kn)
            tjd_uni = np.linspace(tjd[0],tjd[1],Nstep)
            iso_uni = np.array(map(ft,tjd_uni))
            res = [iso_uni]
        elif len(time_kn) == 3:
            tjd = map(lambda x: x.jd, time_kn) 
            tjd1 = np.linspace(tjd[0],tjd[1],Nstep)
            tjd2 = np.linspace(tjd[1],tjd[2],Nstep)
            iso1 = np.array(map(ft,tjd1))
            iso2 = np.array(map(ft,tjd2))
            res = [iso1,iso2] 
        else:
            logging.error('Creation intermediate times needs input length: 2-3')
            exit(1)
        return res


    @classmethod
    def avoid_moon(cls,day_ini,earth_loc):
        #set a radius of N degrees from moon, if full??
        moon1 = apy_coord.get_moon(day_ini,earth_loc)
        return False


    @classmethod
    def point(cls,site_name=None,utc_diff=-3):
        '''This method has few main steps:
        1) gets the site location from another class
        2) gets the range of hours the night lasts
        3) pick the coordinates of one field and gets its AltAz position at
        a given time
        PENDING:
        - if the field is visible, get the airmass (use threshold by user)
        and select the order of observing

        Note:
        - time from user must be feed in local time, but processing is made in 
        UTC. The results are given in LOCAL time by astropy
        '''
        #Earth location
        site = Telescope.site(name=site_name)
        #some SkyCoordinate, from Event 3
        ra,dec = [0.430246,-35.4841]
        #the RA DEC coordinates by convenience must be translated to
        #altitude(0-90) azimuth(0-360,N->E)
        
        #time range: must be set using the Sun elevation
        #now I'll use april 13, from 8.34pm to 6.57am 
        #Chile: UTC-3; Daylight save: UTC-4 -> local = UTC - <3,4>
        #deltaUTC = utc_diff*apy_u.hour
        #local_t1 = apy_time.Time('2017-04-13 20:34:00')
        #local_t2 = apy_time.Time('2017-04-14 06:57:00')
        #trange = [local_t1,local_t2]
        #trange = [x - deltaUTC for x in trange]
        
        #starting at noon, scan for 24hrs searching for the Sun at -14deg 
        deltaUTC = utc_diff*apy_u.hour
        taux = apy_time.Time('2017-04-13 12:00:00') - deltaUTC
        t_window = Schedule.eff_night(taux,site)
        #define times for the range of hours of observation window
        timeshot = Schedule.scan_night(t_window,Nstep=10)
        
        '''if only half nites are used, then timeshot will have 1 array,
        if whole night, 2 arrays
        '''

        #for each of the time stamps, find the zenith RA
        fx = lambda y: Telescope.zenith_ra(y,site)        
        zen_ra = [fx(tm) for tm in timeshot]
        
        #with RA for the zenith and using DEC as it, translate the borders
        #found in CTIO PDF document, given in HourAngle,Dec



        #get the telescope horizon limits, in AltAz coordinates
        #lim_alt,lim_az,lim_altfix = Telescope.horizon_limits_tcs(t_window,site)
        
        #give time range with more entries? yes, because an object can
        #be at the observability window maybe for only a part of the night
        exit()

        #transformation 1: non considering obstime in RADEC initialization
        radec = apy_coord.SkyCoord(ra=ra,dec=dec,frame='icrs',
                                unit=(apy_u.deg,apy_u.deg))
        f = lambda x : radec.transform_to(apy_coord.AltAz(obstime=x,
                                                        location=site))
        #list of objects AltAz
        altaz = [f(a) for a in trange]
        #to transform AltAz to angles in degrees, use Angle(altaz[0].alt).dms

        
        print apy_coord.Angle(altaz[0].alt).dms[:]
        print apy_coord.Angle(altaz[0].alt).is_within_bounds('0d', '360d')
        print altaz[0].alt
        print altaz[0].az
        print altaz[0].secz
        print apy_coord.Angle([-20, 150, 350]*apy_u.deg).is_within_bounds(
            -10*apy_u.deg,None)
        #transformation 2: considering obstime in RADEC
        #Nothing changes, so I will not use the additional arg





if __name__ == '__main__':

    print 'starting...'
    Schedule.point()
    #Telescope.site()
