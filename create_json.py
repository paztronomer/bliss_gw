"""Francisco Paz-Chinchon, University of Illinois/NCSA
DES Collaboration"""
import os
import sys
import time
import subprocess
import shlex
import argparse
import timeit
import math
import numpy as np
import pandas as pd
import astropy
import astropy.units as apy_u
import astropy.coordinates as apy_coord
import astropy.time as apy_time
import scipy.signal as signal #v0.19
import scipy.interpolate as interpolate
import logging
import easyaccess as ea


class Toolbox():
    @classmethod
    def dbquery_ea(cls,expnum_list,band=None,unique_band=False):
        if len(expnum_list) == 1:
            expnum_aux = expnum_list[0]
        elif len(expnum_list) > 1:
            expnum_aux = ",".join(map(str,expnum_list))
        else:
            logging.error("EXPNUM list has zero entries")
            exit(1)
        connect = ea.connect("desoper")
        cursor = connect.cursor()
        if unique_band:
            qaux = " and band='{0}'".format(band)
        else:
            qaux = ""
        q = "select expnum,radeg,decdeg,object,band from exposure"
        q += " where expnum in ({0})".format(expnum_aux)
        q += qaux
        q += " order by expnum"
        #df_obj = cursor.execute(q)
        df_obj = connect.query_to_pandas(q)
        connect.close()
        if (len(df_obj.index) == 0):
            logging.error("\nNo entries in the DB for the inputs EXPNUM")
            exit(1)
        elif (np.unique(df_obj["BAND"].values).shape[0] > 1):
            logging.warning("\nEXPNUMs from different bands")
            print "Bands: {0}".format(np.unique(df_obj["BAND"].values))
        #for unique band, must stopt the code if not matches criteria
        if unique_band and (len(df_obj.index) != len(expnum_list)):
            logging.error("\nMultiple bands for the input EXPNUMs")
            print "Bands: {0}".format(np.unique(df_obj["BAND"].values))
            exit(1)
        return df_obj

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
        Inputs
        - x,y: one dimensional arrays
        - degree: integer, defines the polynomial fitting
        Returns
        - object, the interpolator

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
        t = np.r_[(x[0],)*(degree+1),t,(x[-1],)*(degree+1)]
        lsq_spline = interpolate.make_lsq_spline(x,y,t,degree)
        return lsq_spline

    @classmethod
    def delta_hr(cls,time_arr):
        """Method to return a round amount of hours for a astropy TimeDelta
        object, calculated from the peak-to-peak value given an astropy
        Time array
        Inputs
        - time_arr: array of Time objects
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


class Loader():
    @classmethod
    def obj_field(cls,path,fname):
        """Method to open the tables containing RA,DEC,EXPNUM from the list
        of objects and return a list of pandas tables. Note that inputs
        tables must have a capitalized header with RA, DEC, and EXPNUM on it
        Inputs
        - path: string containing parent path to the tables.
        - fname: list of strings, containing the filenames of the object tables
        Returns
        - list of pandas objects, being the loaded RA DEC tables
        """
        tablst = []
        for fi in fname:
            aux_fn = os.path.join(path,fi)
            tmp = pd.read_table(aux_fn,sep="\s+",usecols=["RA","DEC","EXPNUM"],
                            engine="python",comment="#")
            tablst.append(tmp)
        return tablst


class JSON():
    def __init__(self,
                count=1,
                note="Added to queue by user, not obstac",
                seqid_LIGO=0,seqtot=0,seqnum=0,
                objectname=None,
                propid=None,
                exptype="object",
                progr=None,
                ra=None,dec=None,
                band="i",
                exptime=90,
                til_id=1,
                comment=None,
                towait="False"):
        """Method to fill the dictionary for each of the objects passing the
        observability criteria
        Format:
        - general open/close: []
        - per entry open/close: {}
        - separator: comma"""
        d = dict()
        d["count"] = count
        d["note"] = note
        d["seqid"] = seqid_LIGO #ask as input, use event 1,2,3
        d["seqtot"] = seqtot
        d["seqnum"] = seqnum
        d["object"] = objectname #slack, "DESGW hex (RA)-(DEC) tiling 1"
        d["propid"] = propid #slack
        d["expType"] = exptype
        d["program"] = progr #slack
        d["ra"] = ra
        d["dec"] = dec
        d["filter"] = band
        d["exptime"] = exptime
        d["tiling_id"] = til_id #slack
        d["comment"] = comment
        d["wait"] = towait
        self.dictio = d

    def write_out(self,wrout,itera,maxitera):
        """Method to write the plain text file with the JSON information.
        Receives a writable object and works on it
        """
        if itera == 0:
            wrout.write("[\n")
            wrout.write("\t{\n")
            for idx,(key,value) in enumerate(self.dictio.iteritems()):
                if isinstance(value,str) and np.less(idx,len(self.dictio)-1):
                    wrout.write("\t\t\"{0}\": \"{1}\",\n".format(key,value))
                elif isinstance(value,str) and np.equal(idx,len(self.dictio)-1):
                    wrout.write("\t\t\"{0}\": \"{1}\"\n".format(key,value))
                elif np.less(idx,len(self.dictio)-1):
                    wrout.write("\t\t\"{0}\": {1},\n".format(key,value))
                else:
                    wrout.write("\t\t\"{0}\": {1}\n".format(key,value))
            wrout.write("\t},\n")
        elif itera == maxitera:
            wrout.write("\t{\n")
            for idx,(key,value) in enumerate(self.dictio.iteritems()):
                if isinstance(value,str) and np.less(idx,len(self.dictio)-1):
                    wrout.write("\t\t\"{0}\": \"{1}\",\n".format(key,value))
                elif isinstance(value,str) and np.equal(idx,len(self.dictio)-1):
                    wrout.write("\t\t\"{0}\": \"{1}\"\n".format(key,value))
                elif np.less(idx,len(self.dictio)-1):
                    wrout.write("\t\t\"{0}\": {1},\n".format(key,value))
                else:
                    wrout.write("\t\t\"{0}\": {1}\n".format(key,value))
            wrout.write("\t}\n")
            wrout.write("]\n")
        else:
            wrout.write("\t{\n")
            for idx,(key,value) in enumerate(self.dictio.iteritems()):
                if isinstance(value,str) and np.less(idx,len(self.dictio)-1):
                    wrout.write("\t\t\"{0}\": \"{1}\",\n".format(key,value))
                elif isinstance(value,str) and np.equal(idx,len(self.dictio)-1):
                    wrout.write("\t\t\"{0}\": \"{1}\"\n".format(key,value))
                elif np.less(idx,len(self.dictio)-1):
                    wrout.write("\t\t\"{0}\": {1},\n".format(key,value))
                else:
                    wrout.write("\t\t\"{0}\": {1}\n".format(key,value))
            wrout.write("\t},\n")


class Telescope():
    @classmethod
    def site(cls,name=None):
        """Method to setup a coordinate system for the site of observation.
        By default, CTIO is assumed, given its geodetic coordinates, taken from
        http://www.ctio.noao.edu/noao/node/2085 where:
        - lat/long: [deg,min,sec]
        - altitude: mt
        - ellipsoid: WGS84 where South and West are negative
        If a new site is requested, 'name' must be one of the list in
        EarthLocation... but remember the horizon limits are defined for Blanco
        Telescope
        Inputs
        - name: string representing the name of the observing site, as given in
        astropy.coordinates.EarthLocation.get_site_names(). If none, then uses
        CTIO preloaded data
        Returns
        - coordinate system
        """
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
        (HOUR_ANGLE=0) == (RA at ZENITH)
        Inputs
        - ra_zenith: array of zenith RA, for different times. The times are
        the same as those given in time_interp
        - time_interp: arrays of astropy time objects, over which the
        interpolation is made
        Returns
        - 3 lists: one of the fit, other for time and its RA at zenith, and
        other for DEC limits.

        About the interpolator object, scipy.interpolate._bsplines.BSpline:
        it acts as a function of DEC, f(DEC)=RA. Inside time_ra list, elements
        are: (1) time for this setup, (2) RA of the zenith

        Notes:
        - the interpolator is applied for the positive part of HourAngle
        - When this interpolator is applied, returns an array
        - PLC limits aren't taken into account
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
        dec_lim = [Dec_d.min(),Dec_d.max()]
        #Steps:
        #1) for each of the time stamps, transform HourAngle to RA
        #2) return the interpolator, using RA
        #Note: I will calculate only the upper envelope, as this is symmetric
        fit,time_ra = [],[]
        for idx,zen in enumerate(ra_zenith[:,0]):
            tmp_ra = HAngle.degree + zen
            lsq_dec = Toolbox.lsq_interp(np.array(Dec_d),np.array(tmp_ra))
            fit.append(lsq_dec)
            time_ra.append((time_interp[idx],zen))
        #Usage: xs, tmp_lsq(xs)
        return fit,time_ra,dec_lim

    @classmethod
    def zenith_ra(cls,time_arr,site):
        """Calculates the RA of the azimuth, at a given time, at a given site.
        Time must be given in UTC.
        Inputs
        - array of times (astropy.time class) containing the set of time stamps
        between the borders. The shape is (interpolation_inside_interval,1)
        - site: location of the observing site
        Returns
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

    @classmethod
    def altaz_airm(cls,coord,time,site):
        """For a single position, calculate the airmass at a given time and
        location, given its RA-DEC coordinates
        Inputs
        - coord: tuple containing RA and DEC, both in degrees
        - time: astropy-time object (array), describing the observation
        time
        - site: descriptor for the site (coordinate system)
        Returns
        - a tuple with (altitude,azimuth,airmass)
        """
        radec = apy_coord.SkyCoord(ra=coord[0],dec=coord[1],frame="icrs",
                                unit=(apy_u.deg,apy_u.deg))
        f = lambda x : x.transform_to(apy_coord.AltAz(obstime=time,
                                    location=site))
        res = f(radec)
        return (res.alt.degree[0],res.az.degree[0],res.secz[0])


class Schedule():
    @classmethod
    def eff_night(cls,day_ini,earth_loc):
        """This method calculates the time range between the Sun being
        at -14deg below the horizon (begin and end of night), for a single
        night of observation
        Inputs
        - day_ini: noon of the initial day of observation, the date at which
        the night starts
        - earth_loc: location of the observing site, EarthLocation object
        Returns
        - array of astropy.time.core.Time elements, containing begin,
        middle, and end of the night
        """
        #sample the day by minute
        aux = day_ini + np.linspace(0,24,1440) * apy_u.hour
        #Sun position is in GCRS
        sun1 = apy_coord.get_sun(aux)
        altaz_frame = apy_coord.AltAz(obstime=aux,location=earth_loc)
        sun1_altaz = sun1.transform_to(altaz_frame)
        #transform AltAz to degrees
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
    def scan_night(cls,time_kn,Nstep=100):
        """Method to interpolate the night range, in astropy Time objects.
        Use 2 or 3 times for begin,middle,end of the night to create a
        set of intermediate values. The number of steps given as argument may
        be transparent to the user.
        Inputs
        - time_kn: 1D array with 2 or 3 astropy Time entries, representing
        begin, (optional: middle), and end of the night
        - Nstep: number of steps at which calculate the time interpolation
        Returns
        - array of shape (M,N) where N (N=1) is  the number of sections of the
        night (values can be 2 or 3), and M the number of time intervals the
        section of the night want to be divided for a later interpolation.
        Each time is an object of the class astropy.time.core.Time
        """
        #to iso format
        ft = lambda x: apy_time.Time(apy_time.Time(x,format="jd"),
                                    scale="utc",format="iso")
        #cases for 1 and 2 divisions of the night
        tjd = map(lambda x: x.jd, time_kn)
        if (time_kn.shape[0] in [2,3]):
            tjd_uni = np.linspace(np.min(tjd),np.max(tjd),Nstep+1)
            iso_uni = np.array(map(ft,tjd_uni))
            res = iso_uni
        else:
            logging.error("Scan method needs 2-3 astropy-Time inputs")
            exit(1)
        res = res[:,np.newaxis]
        return res

    @classmethod
    def avoid_moon(cls,day_ini,earth_loc):
        """If implemented, use apy_coord.get_moon(day_ini,earth_loc)
        """
        pass

    @classmethod
    def point_onenight(cls,
                    path_tab=None,
                    fname_csv=None,
                    fname_json=None,
                    object_list=None,
                    site_name=None,
                    utc_minus_local=None,
                    begin_day=None,
                    obs_interval=None,
                    T_step=None,
                    max_airm=None,
                    count=None,seqid_LIGO=None,propid=None,exptype=None,
                    progr=None,band=None,exptime=None,til_id=None,
                    comment=None,note=None,towait=None,unique_band=None):
        """Method to wrap different methods, to calculate objects observability
        for a single night

        This method has few main steps:
        1) gets the site location from the Telescope class
        2) gets the range of hours every night lasts (to get the splitting for
        half observation nights). Remember the transformation to UTC so all
        results are not in local time.
        3) inside the observation window, define a grid of times to be used as
        steps for the calculations of objects observability
        4) given the borders or the observation window for a specific night,
        interpolate as many points as the interpolation rate was defined. This
        is made inside the interval, as defined by the input borders.
        5) calculate the RA of the zenith in each of the time stamps
        6) giving RA for the zenith at different times, locate the CTIO horizon
        limits. Translate the borders found in CTIO PDF document to RA-DEC,
        initially given in HourAngle,Dec
        7) for each GW event, for each object coordinate, for each of the
        time stamps (inside one night), see if fit inside bounds. Then, if
        inside bounds, see if the position matches the ALtAz airmass criteria
        8) within the selection, we can have multiple entries per night
        (different time stamps of each night subsampling), so must decide
        based in the lowest angle from zenith
        9) after pick the lowest airmass, we do have the results for the object
        matching the criteria. Then write out both the results and the JSON
        files. A different JSON file for each date

        Note:
        - time from user must be feed in local time, but processing is made in
        UTC. The results are given in UTC time by astropy
        """
        site = Telescope.site(name=site_name)
        #starting at noon, scan for 24hrs searching for the Sun at -14deg,
        #returns an array on which time entries are astropy.time.core.Time
        deltaUTC = utc_minus_local*apy_u.hour
        delta_noon = 12.0*apy_u.hour
        t_utc = apy_time.Time(begin_day) + deltaUTC + delta_noon
        t_window = Schedule.eff_night(t_utc,site)
        if obs_interval == "first":
            t_window = t_window[:-1]
        elif obs_interval == "second":
            t_window = t_window[1:]
        elif (obs_interval == "full") or (obs_interval is None):
            pass
        else:
            logging.error("Interval must be: first, second, full")
        #Returns an array of astropy.time.core.Time entries
        N_hr,N_min = Toolbox.delta_hr(t_window)
        if T_step is None:
            xstep = np.int(np.round(N_min/10.))
        elif np.less(T_step,N_min):
            xstep = np.int(np.round(N_min/np.float(T_step)))
        else:
            logging.error("Time step (T_step) must be smaller than the window")
            exit(1)
        t_interp = Schedule.scan_night(t_window,Nstep=xstep)

        #Returns an array of same shape as the array of interpolated times
        zen_ra = Telescope.zenith_ra(t_interp,site)

        #output: [object f(dec)=ra],[time,RA of zenith],[min(dec),max(dec)]
        func_dec,time_RA,declim = Telescope.horizon_limits_tcs(zen_ra,t_interp)

        #iterate over list of objects
        radec = Loader.obj_field(path_tab,object_list)
        sel = []
        for df in radec:
            #get the object names from desoper DB
            dbinfo = Toolbox.dbquery_ea(list(df["EXPNUM"].values),band=band,
                                    unique_band=unique_band)
            for index,row in df.iterrows():
                for idx0,tRA in enumerate(time_RA):
                    low_RA = tRA[1] - func_dec[idx0](row["DEC"])
                    cond1 = np.less_equal(row["RA"],func_dec[idx0](row["DEC"]))
                    cond2 = np.greater_equal(row["RA"],low_RA)
                    cond3 = np.less_equal(row["DEC"],declim[1])
                    cond4 = np.greater_equal(row["DEC"],declim[0])
                    if cond1 and cond2 and cond3 and cond4:
                        args = [(row["RA"],row["DEC"]),tRA[0],site]
                        alt,az,secz = Telescope.altaz_airm(*args)
                        cond5 = np.less_equal(secz,max_airm)
                        cond6 = np.greater_equal(secz,1.)
                        if cond5 and cond6:
                            dfaux = dbinfo.loc[dbinfo["EXPNUM"]==row["EXPNUM"]]
                            z_tmp = math.degrees(math.acos(1/np.float(secz)))
                            tmp = (index+1,row["RA"],row["DEC"],alt,az)
                            tmp += (z_tmp,secz,tRA[0][0]-deltaUTC)
                            tmp += (dfaux["OBJECT"].values[0],)
                            sel.append(tmp)
        #if no object meets observability criteria
        if len(sel) == 0:
            err_mssg = "\tNo object matches the observability criteria!"
            err_mssg += "\n\tThere will be no output file\n"
            logging.error(err_mssg)
        else:
            #create a pandas DataFrame or structured array for easier selection
            cols = ["n_id","ra","dec","alt","az","z","secz","local_time","obj"]
            sel_df = pd.DataFrame(sel,columns=cols)
            min_df = pd.DataFrame()
            for N in sel_df["n_id"].unique():
                tmp_df = sel_df.loc[(sel_df["n_id"]==N)]
                tmp_df = tmp_df.loc[tmp_df["secz"].idxmin(axis="columns")]
                min_df = min_df.append(tmp_df)
            #sort by RA and then by DEC
            min_df.sort_values("secz",ascending=True,inplace=True)
            #reset index to avoid confusion
            min_df = min_df.reset_index()
            #write out the resume table
            min_df.to_csv(fname_csv,sep=",",index=False,header=True)
            #write json file for this night
            fjson = open(fname_json,"w")
            for index,row in min_df.iterrows():
                jw = dict()
                jw["count"] = count
                jw["seqid_LIGO"] = seqid_LIGO
                jw["propid"] = propid
                jw["exptype"] = exptype
                jw["progr"] = progr
                jw["band"] = band
                jw["exptime"] = exptime
                jw["til_id"] = til_id
                jw["comment"] = comment
                jw["note"] = note
                jw["towait"] = towait
                jw["seqtot"] = len(min_df.index)
                jw["seqnum"] = index + 1
                jw["objectname"] = row["obj"]
                jw["ra"] = row["ra"]
                jw["dec"] = row["dec"]
                JSON(**jw).write_out(fjson,index,len(min_df.index)-1)
            fjson.close()

    @classmethod
    def point_allnight(cls,
                    path_tab=None,
                    path_out=None,
                    root_csv=None,
                    root_json=None,
                    date_tab=None,
                    object_list=None,
                    utc_minus_local=4,
                    T_step=20,
                    max_airm=1.8,
                    count=None,seqid_LIGO=None,propid=None,exptype=None,
                    progr=None,band=None,exptime=None,til_id=None,
                    comment=None,note=None,towait=None,unique_band=None):
        """Method to iteratively call the observability calculation, night by
        night
        """
        if root_csv is None:
            root_csv = "sel_info"
        if root_json is None:
            root_json = "obs"
        if propid is None:
            logging.warning("\t(!) No value has been entered for Proposal ID")
        if unique_band is None:
            unique_band = False
        #
        date_fn = os.path.join(path_tab,date_tab)
        wd = pd.read_table(date_fn,sep="\s+",names=["date","part"],
                        header=None,engine="python",comment="#")
        #
        for idx,row in wd.iterrows():
            print "Working on night: {0}/{1} ".format(*row)
            t0 = time.time()
            out_aux = row["date"][5:].replace("-","_")
            #csv resume table filename
            out_csv = "{0}_{1}.csv".format(root_csv,out_aux)
            out_csv = os.path.join(path_out,out_csv)
            #json file filename
            out_json = "{0}_{1}.json".format(root_json,out_aux)
            out_json = os.path.join(path_out,out_json)
            #arguments for night to night call
            kw = dict()
            kw["path_tab"] = path_tab
            kw["fname_csv"] = out_csv
            kw["fname_json"] = out_json
            kw["object_list"] = object_list
            kw["utc_minus_local"] = utc_minus_local
            kw["begin_day"] = row["date"]
            kw["obs_interval"] = row["part"]
            kw["T_step"] = T_step
            kw["max_airm"] = max_airm
            kw["count"] = count
            kw["seqid_LIGO"] = seqid_LIGO
            kw["propid"] = propid
            kw["exptype"] = exptype
            kw["progr"] = progr
            kw["band"] = band
            kw["exptime"] = exptime
            kw["til_id"] = til_id
            kw["comment"] = comment
            kw["note"] = note
            kw["towait"] = towait
            kw["unique_band"] = unique_band
            Schedule.point_onenight(**kw)
            t1 = time.time()
            print "Elapsed time: {0:.2f} minutes".format((t1-t0)/60.)


if __name__ == "__main__":
    print "\tRunning script: {0}\n\t{1}".format(__file__,"="*30)

    """Fill up the different arguments from the command line call
    """
    gral_descr = "Script to calculate observability from CTIO, using"
    gral_descr += " Blanco 4m telecope setup, for a set (or single) of "
    gral_descr += " coordinates given the observing date. There are 2 type of"
    gral_descr += " optional arguments: those wo serve as setup and those who"
    gral_descr += " will be directly inserted in the JSON files. NOTE: for"
    gral_descr += " JSON arguments use quotes if any space is in the string"
    aft = argparse.ArgumentParser(description=gral_descr)
    #positional
    h1 = "Set of space-separated filenames for the tables containig the"
    h1 += " coordinates to be calculated. Format: tables must have RA, DEC, "
    h1 += " and EXPNUM in header, despite the other column names; separator"
    h1 += " is expected to be spaces"
    aft.add_argument("objects",help=h1,nargs="+")
    h2 = "Unique file containig all the nights for which the observability"
    h2 += " will be calculated. Format: 2 columns being YYYY-MM-DD"
    h2 += " {first/second/all}"
    aft.add_argument("dates",help=h2)
    #optional
    h3 = "Path to folder containing the source tables for objects and dates."
    h3 += " Default is current directory"
    aft.add_argument("--source","-s",help=h3,metavar="",default=os.getcwd())
    h4 = "Path to folder for output files (JSON and CSV)."
    h4 += " Default is current directory"
    aft.add_argument("--out","-o",help=h4,metavar="",default=os.getcwd())
    h5 = "Root string of output name(s) for the file(s) containing"
    h5 += " the coordinates that passed the observability criteria, plus"
    h5 += " additional information. Default is \'sel_info\'"
    aft.add_argument("--root_csv",help=h5,metavar="",default="selected")
    h6 = "Root string of JSON output files (if objects were found)."
    h6 += " Default is \'obs\'"
    aft.add_argument("--root_json",help=h6,metavar="",default="obs")
    h7 = "Difference in hours between UTC and the observing location, namely,"
    h7 += " CTIO (value=UTC-LOCAL). Default: 4"
    aft.add_argument("--utc_diff","-u",help=h7,type=float,metavar="",default=4)
    h8 = "Step (in minutes) at which the night will be sampled to calculate"
    h8 += " the observability at each interval. Default: 10"
    aft.add_argument("--t_step","-t",help=h8,metavar="",type=float,default=10)
    h9 = "Maximum airmass at which the objects want to be observed."
    h9 += " Default: 1.8"
    aft.add_argument("--max_airmass","-m",help=h9,metavar="",default=1.8,
                    type=float)
    h21 = "Whether to request all input EXPNUMs to have the same band as the"
    h21 += " given in '--band' argument (where deafult is 'i')"
    aft.add_argument("--req_one_band",help=h21,action='store_true')
    #optional to be added in JSON files
    h10 = "JSON optional.Number of exposures to be taken for each object."
    h10 += " Default: 1"
    aft.add_argument("--count",help=h10,metavar="",default=1,type=int)
    h11 = "JSON required. Sequence ID, eg: \'LIGO event x\'"
    aft.add_argument("--sequence",help=h11,metavar="")
    h12 = "JSON required. Proposal ID"
    aft.add_argument("--proposal",help=h12,metavar="")
    h13 = "JSON optional. Exposure type. Default: \'object\'"
    aft.add_argument("--exptype",help=h13,metavar="",default="object")
    h14 = "JSON required. Program ID, example: \'BLISS\'"
    aft.add_argument("--program",help=h14,metavar="")
    h15 = "JSON optional. Band to be used. Default: i"
    aft.add_argument("--band",help=h15,metavar="",default="i")
    h16 = "JSON optional. Exposure time in seconds. Default: 90"
    aft.add_argument("--exptime",help=h16,metavar="",default=90,type=float)
    h17 = "JSON optional. ID of the tiling. Default: 1"
    aft.add_argument("--tiling",help=h17,metavar="",default=1)
    h18 = "JSON optional. Comment to be added. Default: empty"
    aft.add_argument("--comment",help="JSON Optional",metavar="",default="")
    h19 = "JSON optional. Note to be added. Default: \'Added to queue by"
    h19 += " user, not obstac\'"
    aft.add_argument("--note",help=h19,metavar="",
                    default="Added to queue by user, not obstac")
    h20 = "JSON optional. Whether to wait to proceed for next exposure."
    h20 += " Default: \'False\'"
    aft.add_argument("--wait",help=h20,metavar="",default="False")
    #parser
    args = aft.parse_args()
    kw0 = vars(args)
    kw1 = dict()
    kw1["object_list"] = kw0["objects"]
    kw1["date_tab"] = kw0["dates"]
    kw1["path_tab"] = kw0["source"]
    kw1["path_out"] = kw0["out"]
    kw1["root_csv"] = kw0["root_csv"]
    kw1["root_json"] = kw0["root_json"]
    kw1["utc_minus_local"] = kw0["utc_diff"]
    kw1["T_step"] = kw0["t_step"]
    kw1["max_airm"] = kw0["max_airmass"]
    kw1["count"] = kw0["count"]
    kw1["seqid_LIGO"] = kw0["sequence"]
    kw1["propid"] = kw0["proposal"]
    kw1["exptype"] = kw0["exptype"]
    kw1["progr"] = kw0["program"]
    kw1["band"] = kw0["band"]
    kw1["exptime"] = kw0["exptime"]
    kw1["til_id"] = kw0["tiling"]
    kw1["comment"] = kw0["comment"]
    kw1["note"] = kw0["note"]
    kw1["towait"] = kw0["wait"]
    kw1["unique_band"] = kw0["req_one_band"]
    #
    Schedule.point_allnight(**kw1)

    if False:
        stp = "from __main__ import Schedule; sch = Schedule"
        timeval = timeit.timeit("sch.point()",setup=stp,number=50)
        print timeval
