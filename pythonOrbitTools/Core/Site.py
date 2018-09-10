##
# @file Site.py
# @brief 
# @author df_justforfun@163.com
# @version 1.0
# @date 2018-07-22

import math
from pythonOrbitTools.Core.Globals import Globals
from pythonOrbitTools.Core.Coord import Geo, TopoTime
from pythonOrbitTools.Core.Eci import EciTime
from pythonOrbitTools.Core.Julian import Julian
from pythonOrbitTools.Core.Vector import Vector

##
# @brief The Site class encapsulates a location on earth.
class Site(object):

    # region Properties

    @property
    ##
    # @brief The name of the location
    #
    # @return 
    def Name(self):
        return self._name

    @property
    ##
    # @brief Latitude, in radians. A negative value indicates latitude south.
    #
    # @return 
    def LatitudeRad(self):
        return self._geo.LatitudeRad

    @property
    ##
    # @brief Longitude, in radians. A negative value indicates longitude west.
    #
    # @return 
    def LongitudeRad(self):
        return self._geo.LongitudeRad

    @property
    ##
    # @brief Latitude, in degrees. A negative value indicates latitude south.
    #
    # @return 
    def LatitudeDeg(self):
        return self._geo.LatitudeDeg

    @property
    ##
    # @brief Longitude, in degrees. A negative value indicates longitude west.
    #
    # @return 
    def LongitudeDeg(self):
        return self._geo.LongitudeDeg

    @property
    ##
    # @brief The altitude of the site above the ellipsoid model, in kilometers.
    #
    # @return 
    def Altitude(self):
        return self._geo.Altitude

    @property
    ##
    # @brief The contained geodetic coordinates.
    #
    # @return 
    def Geo(self):
        return self._geo

    # endregion

    ##
    # @brief Creates a new instance of the class None Geo and None Name
    #
    # @return 
    def __init__(self):
        self._geo  = None
        self._name = None

    ##
    # @brief Initialize the instance of the class with the given components.
    #
    # @param degLat Latitude in degrees (negative south).
    # @param degLon Longitude in degrees (negative west).
    # @param kmAlt Altitude in kilometers.
    # @param name The name of the location, default value is "".
    #
    # @return 
    def InitializeByDegLatAndDegLonAndKmAltAndName(self, degLat, degLon, kmAlt, name=""):
        self._geo = Geo().InitializeByRadLatAndRadLonAndKmAlt(Globals.ToRadians(degLat),
                                                              Globals.ToRadians(degLon),
                                                              kmAlt)
        self._name = name
        return self

    ##
    # @brief Initialize the instance of the class with the given Geo object. 
    #
    # @param geo The Geo object.
    #
    # @return 
    def InitializeByGeo(self, geo):
        self._geo = Geo().InitializeByGeo(geo)
        return self

    ##
    # @brief Calculates the ECI coordinates of the site
    #
    # @param date Time of position calculation.
    #
    # @return The site`s ECI coordinates at the given time.
    def PositionEciByJulianTime(self,date):
        return EciTime().InitializeByGeoAndDate(self._geo, date)

    ##
    # @brief Calculates the ECI coordinates of the site.
    #
    # @param utc Time of position calculation.
    #
    # @return The site`s ECI coordinates at the given time. 
    def PositionEciByDateTime(self, utc):
        return EciTime().InitializeByGeoAndDate(self._geo, Julian().InitializeByUTC(utc))


    ##
    # @brief Returns the topo-centric (azimuth, elevation, etc.) coordinates for
    #        a target object described by the given ECI coordinates.
    #
    # @param eci The ECI coordinates of the target object.
    #
    # @return The look angle to the target object.
    def GetLookAngle(self, eci):
        # Calculate the ECI coordinates for this Site object at the time of interest
        date = eci.Date
        eciSite = self.PositionEciByJulianTime(date)
        vecRgRate = Vector(eci.Velocity.X - eciSite.Velocity.X,
                           eci.Velocity.Y - eciSite.Velocity.Y,
                           eci.Velocity.Z - eciSite.Velocity.Z)

        x = eci.Position.X - eciSite.Position.X
        y = eci.Position.Y - eciSite.Position.Y
        z = eci.Position.Z - eciSite.Position.Z
        w = math.sqrt(Globals.Sqr(x) + Globals.Sqr(y) + Globals.Sqr(z))

        vecRange = Vector(x, y, z, w)

        # The site`s Local Mean Sidereal Time at the time of interest.
        theta = date.ToLmst(self.LongitudeRad)

        sin_lat     = math.sin(self.LatitudeRad)
        cos_lat     = math.cos(self.LatitudeRad)
        sin_theta   = math.sin(theta)
        cos_theta   = math.cos(theta)

        top_s = sin_lat * cos_theta * vecRange.X +\
                sin_lat * sin_theta * vecRange.Y -\
                cos_lat * vecRange.Z
        top_e = -sin_theta * vecRange.X + cos_theta * vecRange.Y
        top_z = cos_lat * cos_theta * vecRange.X +\
                cos_lat * sin_theta * vecRange.Y +\
                sin_lat * vecRange.Z
        az    = math.atan(-top_e/top_s)

        if top_s > 0.0:
            az += Globals.Pi

        if az < 0.0:
            az += 2.0 * Globals.Pi

        el      = math.asin(top_z/vecRange.W)
        rate    = (vecRange.X * vecRgRate.X +\
                   vecRange.Y * vecRgRate.Y +\
                   vecRange.Z * vecRgRate.Z)/vecRange.W

        topo = TopoTime().InitializeByRadAzAndRadElAndRangeAndRangeRateAndDate(az,
                                                                               el,
                                                                               vecRange.W,
                                                                               rate,
                                                                               eci.Date)

        # if WANT_ATMOSPHERIC_CORRECTION
        # endif

        return topo



    ##
    # @brief Converts to a string representation of the form "120.00N 90.00W 500m"
    #
    # @return The formatted string.
    def __str__(self):
        return self._geo.__str__()
