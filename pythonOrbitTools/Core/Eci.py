##
# @file Eci.py
# @brief 
# @author df_justforfun@163.com
# @version 1.0
# @date 2018-07-19

import math
from pythonOrbitTools.Core.Globals import Globals
from pythonOrbitTools.Core.Vector import Vector

##
# @brief Encapsulates an Earth-Centered Inertial coordinate position/velocity
class Eci(object):

    # region Properties

    @property
    def Position(self):
        return self._pos

    @property
    def Velocity(self):
        return self._vel

    # endregion

    ##
    # @brief: Creates a new instance of the class None position and None velocity. 
    #
    # @return 
    def __init__(self):
        self._pos = None
        self._vel = None

    ##
    # @brief: Initialize the instance of the class with the given position and 
    #         velocity components.
    #
    # @param pos: The position vector.
    # @param vel: The velocity vector.
    #
    # @return 
    def InitializeByPosAndVel(self, pos, vel):
        self._pos = pos
        self._vel = vel

        return self

    ##
    # @brief: Initialize the isinstance of the class with the given geodetic coordinates. 
    #         Assumes the Earth is an oblate spheroid.
    #         Reference: The 1992 Astronomical Almanac, page K11
    #         Reference: www.celestrak.com (Dr. T.S. Kelso)
    #
    # @param geo: The geocentric coordinates.
    # @param date: The Julian date.
    #
    # @return 
    def InitializeByGeoAndDate(self, geo, date):
        lat = geo.LatitudeRad
        lon = geo.LongitudeRad
        alt = geo.Altitude

        # Calculate Local Mean Sidereal Time (theta)
        theta = date.ToLmst(lon)
        c     = 1.0 / math.sqrt(1.0+Globals.F*(Globals.F-2.0)*Globals.Sqr(math.sin(lat)))
        s     = Globals.Sqr(1.0-Globals.F)*c
        achcp = (Globals.Xkmper*c+alt)*math.cos(lat)

        x = achcp * math.cos(theta) # km
        y = achcp * math.sin(theta) # km
        z = (Globals.Xkmper*s+alt)*math.sin(lat) # km
        w = math.sqrt(Globals.Sqr(x)+Globals.Sqr(y)+Globals.Sqr(z)) # range, km
        self._pos = Vector(x, y, z, w)

        mfactor = Globals.TwoPi*(Globals.OmegaE/Globals.SecPerDay)
        x = -mfactor*self._pos.Y # km / sec 
        y =  mfactor*self._pos.X # km / sec
        z = 0.0                  # km / sec
        w = math.sqrt(Globals.Sqr(x)+Globals.Sqr(y)) # range rate km/sec^2
        self._vel = Vector(x, y, z, w)


    ##
    # @brief: Scale the position vector by a factor.
    #
    # @param factor
    #
    # @return 
    def ScalePosVector(self, factor):
        self._pos.Mul(factor)

    ##
    # @brief: Scale the velocity vector by a factor. 
    #
    # @param factor
    #
    # @return 
    def ScaleVelVector(self, factor):
        self._vel.Mul(factor)


    ##
    # @brief Returns a string representation of the coordinate and velocity XYZ values.
    #
    # @return: The formatted string. 
    def __str__(self):
        return "km: ({}, {}, {}) km/s:({}, {}, {})".format(self._pos.X,
                                                           self._pos.Y,
                                                           self._pos.Z,
                                                           self._vel.X,
                                                           self._vel.Y,
                                                           self._vel.Z)

##
# @brief: Encapsulates an Earth Centered Inertial coordinate and associated time. 
class EciTime(Eci):

    # region Properties

    @property
    ##
    # @brief The Julian date associated with the ECI coordinates.
    #
    # @return 
    def Date(self):
        return self._date

    # endregion

    ##
    # @brief: Creates an instance of the class None date 
    #
    # @return 
    def __init__(self):
        super().__init__()
        self._date = None


    ##
    # @brief: Initialize the instance of the class with the given position, velocity, and Julian date. 
    #
    # @param pos: The position vector.
    # @param vel: The velocity vector.
    # @param date: The Julian date associated with the position.
    #
    # @return 
    def InitializeByPosAndVelAndDate(self, pos, vel, date):
        super().InitializeByPosAndVel(pos, vel)
        self._date = date
        return self

    ##
    # @brief: Initialize the instance of the class with the given ECI-time coordinates. 
    #
    # @param eci: The ECI coordinates.
    # @param date: The Julian date associated with the ECI coordinates.
    #
    # @return 
    def InitializeByEciAndDate(self, eci, date):
        self.InitializeByPosAndVelAndDate(eci.Position, eci.Velocity, date)
        return self

    ##
    # @brief: Initialize the instance of the class with the given geodetic coordinates and Julian date. 
    #
    # @param geo: The geodetic coordinates.
    # @param date: The Julian date associated with the ECI coordinates.
    #
    # @return 
    def InitializeByGeoAndDate(self, geo, date):
        super().InitializeByGeoAndDate(geo, date)
        self._date = date
        return self

    ##
    # @brief: Initialize the instance of the class with the given geodetic-time coordinates. 
    #
    # @param geotime: The geodetic-time coordinates.
    #
    # @return 
    def InitializeByGeoTime(self, geotime):
        self.InitializeByGeoAndDate(geotime, geotime.Date)
        return self

if __name__ == "__main__":
    pass

