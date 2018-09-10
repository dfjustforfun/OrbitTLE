
##
# @file Globals.py
# @brief Numerical constants
# @author df_justforfun@163.com
# @version 1.0
# @date 2018-07-17

import math

class Globals(object):

    # region Constants

    Pi            = 3.141592653589793
    TwoPi         = 2.0 * Pi
    RadsPerDegree = Pi / 180.0
    DegreesPerRad = 180.0 / Pi

    Gm            = 398601.2 # Earth gravitational constant, km^3/sec^2
    GeoSyncAlt    = 42241.892 # km
    EarthDiam     = 12800.0 # km
    DaySidereal   = (23*3600)+(56*60)+4.09 # sec
    DaySolar      = (24*3600) # sec

    Ae            = 1.0
    Au            = 149597870.0 # Astronomical unit (km) (IAU 76)
    Sr            = 696000.0 # Solar radius (km) (IAU 76)
    Xkmper        = 6378.135 # Earth equatorial radius - kilometers (WGS '72)
    F             = 1.0 / 298.26 # Earth flattening (WGS '72)
    Ge            = 398600.8 # Earth gravitational constant (WGS '72)
    J2            = 1.0826158E-3 # J2 harmonic (WGS '72)
    J3            = -2.53881E-6 # J3 harmonic (WGS '72)
    J4            = -1.65597E-6 # J4 harmonic (WGS '72)
    Ck2           = J2 / 2.0
    Ck4           = -3.0 * J4 / 8.0
    Xj3           = J3
    Qo            = Ae + 120.0 / Xkmper # Globals.Xkmper
    S             = Ae + 78.0 / Xkmper # Globals.Xkmper
    HoursPerDay   = 24.0
    MinPerDay     = 1440.0 # Minutes per day (solar)
    SecPerDay     = 86400.0 # Seconds per day (solar)
    OmegaE        = 1.00273790934 # Earth rotation per sidereal day
    Xke           = math.sqrt(3600.0 * Ge / (Xkmper*Xkmper*Xkmper)) # sqrt(ge) ER^3/min^2
    Qoms2t        = math.pow((Qo-S), 4) # (Qo-S)^4 ER^4

    # endregion

    # region Utility

    @staticmethod
    def Sqr(x):
        return x * x

    @classmethod
    def Fmod2p(cls, arg):
        modu = (arg % cls.TwoPi)

        if modu < 0.0:
            modu += cls.TwoPi

        return modu

    @classmethod
    ##
    # @brief ArcTangent of sin(x) / cos(x). The advantage of this function over arctan()
    #        is that it returns the correct quadrant of the angle.
    #
    # @param cls
    # @param sinx
    # @param cosx
    #
    # @return 
    def AcTan(cls, sinx, cosx):

        if cosx == 0.0:
            if sinx > 0.0:
                ret = cls.Pi / 2.0
            else:
                ret = 3.0 * cls.Pi / 2.0
        else:
            if cosx > 0.0:
                ret = math.atan(sinx/cosx)
            else:
                ret = cls.Pi + math.atan(sinx/cosx)

        return ret

    @classmethod
    def ToDegrees(cls, radians):
        return radians * cls.DegreesPerRad

    @classmethod
    def ToRadians(cls, degrees):
        return degrees * cls.RadsPerDegree

    # endregion

if __name__ == "__main__":

    print("******************* Globals.Sqr **************************")
    print("Globals.Sqr({}) = {}\n".format(3, Globals.Sqr(3)))

    print("******************* Globals.ToDegrees **************************")
    print("Globals.ToDegrees({}) = {} deg\n".format("Pi rad", Globals.ToDegrees(Globals.Pi)))

    print("******************* Globals.ToRadians **************************")
    print("Globals.ToRadians({}) = {} rad\n".format("180 deg", Globals.ToRadians(180)))

    print("******************* Globals.Ck2 **************************")
    print("Globals.Ck2 = {}\n".format(Globals.Ck2))
