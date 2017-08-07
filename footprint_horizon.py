"""Code to plot the footprint polygon plus CTIO horizon limits
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
from astropy import units as u
from astropy.coordinates import SkyCoord
import pandas as pd


class Polygon():
    @classmethod
    def poly(cls):
        #the coordinates for gw1,gw2,gw3
        ra_x,dec_x = [],[]
        for fnm in ["event1_ccds.txt","event2.csv","event3modif_ccds.txt"]:
            df = pd.read_table(fnm,sep="\s+",usecols=["RA","DEC","EXPNUM"],
                            engine="python",comment="#")
            ra_x += list(df["RA"].values)
            dec_x += list(df["DEC"].values)
        c0 = SkyCoord(ra=ra_x*u.degree,dec=dec_x*u.degree,frame='icrs')
        ra_x_rad = c0.ra.wrap_at(180 * u.deg).radian
        dec_x_rad = c0.dec.radian

        #plot the polygon of the footprint
        ra_p,dec_p = np.loadtxt("round13-poly.txt",unpack=True)
        c1 = SkyCoord(ra=ra_p*u.degree,dec=dec_p*u.degree,frame='icrs')
        ra_p_rad = c1.ra.wrap_at(180 * u.deg).radian
        dec_p_rad = c1.dec.radian

        #the horizon
        aux_ra = 48.
        dec_h = np.linspace(-89,37,100)
        # Results of RA is in hours, must be changed to degrees
        ra_h1 = aux_ra + Polygon.horizon()(dec_h)
        ra_h2 = aux_ra - Polygon.horizon()(dec_h)
        c2 = SkyCoord(ra=ra_h1*u.h,dec=dec_h*u.degree,frame='icrs')
        c3 = SkyCoord(ra=ra_h2*u.h,dec=dec_h*u.degree,frame='icrs')
        ra_h1_rad = c2.ra.wrap_at(180 * u.deg).radian
        ra_h2_rad = c3.ra.wrap_at(180 * u.deg).radian
        dec_h_rad = c2.dec.radian


        fig = plt.figure(figsize=(8,5))
        ax1 = fig.add_subplot(111,projection="aitoff")

        plt.grid(True)
        plt.plot(ra_p_rad,dec_p_rad,'-', alpha=0.9)
        plt.plot(ra_h1_rad,dec_h_rad,"g-")
        plt.plot(ra_h2_rad,dec_h_rad,"g-")
        plt.plot(ra_x_rad,dec_x_rad,"ro",markersize=3)

        plt.savefig("fp_horizon.pdf",dpi=150,facecolor="w",edgecolor="w",
                orientation="portrait", papertype=None, format="pdf",
                transparent=False, bbox_inches=None, pad_inches=0.1,
                frameon=None)
        plt.show()

    @classmethod
    def horizon(cls):
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
        #Steps:
        #1) for each of the time stamps, transform HourAngle to RA
        #2) return the interpolator, using RA
        #Note: I will calculate only the upper envelope, as this is symmetric
        tmp_ra = np.array(HAngle)
        lsq_dec = Polygon.lsq_interp(np.array(Dec_d),tmp_ra)
        #usage f(dec)
        return lsq_dec

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

if __name__ == "__main__":
    Polygon.poly()
