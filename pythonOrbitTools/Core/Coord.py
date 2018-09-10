##
# @file Coord.py
# @brief 
# @author df_justforfun@163.com
# @version 1.0
# @date 2018-07-20

import math
from pythonOrbitTools.Core.Globals import Globals

##
# @brief Class to encapsulate geocentric coordinates.
class Geo(object):

    # region Properties

    @property
    ##
    # @brief Latitude, in radians. A negative value indicates latitude south.
    #
    # @return 
    def LatitudeRad(self):
        return self._latitudeRad

    @property
    ##
    # @brief Longitude, in radians. A negative value indicates longitude west.
    #
    # @return 
    def LongitudeRad(self):
        return self._longitudeRad

    @property
    ##
    # @brief Latitude, in degrees. A negative value indicates latitude south.
    #
    # @return 
    def LatitudeDeg(self):
        return Globals.ToDegrees(self._latitudeRad)

    @property
    ##
    # @brief Longitude, in degrees. A negative value indicates longitude west.
    #
    # @return 
    def LongitudeDeg(self):
        return Globals.ToDegrees(self._longitudeRad)

    @property
    ##
    # @brief Altitude, in kilometers, above the ellipsoid model.
    #
    # @return 
    def Altitude(self):
        return self._altitude


    # endregion

    ##
    # @brief Creates an instance of the class None LatitudeRad, None LongitudeRad and None Altitude
    #
    # @return 
    def __init__(self):
        self._latitudeRad  = None
        self._longitudeRad = None
        self._altitude     = None

    ##
    # @brief Initialize the instance of the class with the given source Geo object.
    #
    # @param geo The source Geo object.
    #
    # @return 
    def InitializeByGeo(self, geo):
        self._latitudeRad  = geo.LatitudeRad
        self._longitudeRad = geo.LongitudeRad
        self._altitude     = geo.Altitude
        return self

    ##
    # @brief Initialize the instance of the class with the given components. 
    #
    # @param radLat Latitude, in radians. Negative value indicate latitude south.
    # @param radLon Longitude, in radians. Negative value indicate longitude west.
    # @param kmAlt Altitude above the ellipsoid model, in kilometers.
    #
    # @return 
    def InitializeByRadLatAndRadLonAndKmAlt(self, radLat, radLon, kmAlt):
        self._latitudeRad  = radLat
        self._longitudeRad = radLon
        self._altitude     = kmAlt
        return self

    ##
    # @brief Initialize the instance of the class with the given ECI coordinates.
    #
    # @param eci The ECI coordinates.
    # @param date The Julian date.
    #
    # @return 
    def InitializeByEciAndDate(self, eci, date):
        self.InitializeByPosAndTheta(eci.Position, (Globals.AcTan(eci.Position.Y, eci.Position.X) - date.ToGmst()) % Globals.TwoPi)
        return self


    ##
    # @brief Initialize the instance of the class with the given XYZ coordinates.
    #
    # @param pos
    # @param theta
    #
    # @return 
    def InitializeByPosAndTheta(self, pos, theta):

        theta = theta % Globals.TwoPi

        if theta < 0.0:
            # "wrap" negative modulo
            theta += Globals.TwoPi

        r   = math.sqrt(Globals.Sqr(pos.X) + Globals.Sqr(pos.Y))
        e2  = Globals.F * (2.0 - Globals.F)
        lat = Globals.AcTan(pos.Z, r)

        DELTA = 1.0e-07

        while True:

            phi = lat
            c   = 1.0 / math.sqrt(1.0 - e2 * Globals.Sqr(math.sin(phi)))
            lat = Globals.AcTan(pos.Z + Globals.Xkmper * c * e2 * math.sin(phi), r)

            if math.fabs(lat - phi) <= DELTA:
                break

        self._latitudeRad  = lat
        self._longitudeRad = theta
        self._altitude     = (r/math.cos(lat)) - Globals.Xkmper * c

        return self



    ##
    # @brief Converts to a string representation of the form "38.0N 45.0W 500m"
    #
    # @return 
    def __str__(self):
        latNorth = (self._latitudeRad  >= 0.0)
        lonEast  = (self._longitudeRad >= 0.0)

        # latitude in degrees
        string  = "{}{} ".format(math.fabs(self._latitudeRad), "N" if latNorth else "S")

        # longitude in degrees
        string += "{}{} ".format(math.fabs(self._longitudeRad), "E" if lonEast else "W")

        # elevation in meters
        string += "{}{}m".format(self._altitude*1000.0)

        return string


##
# @brief Class to encapsulate a geocentric coordinate and associated time.
class GeoTime(Geo):

    # region Properties

    @property
    def Date(self):
        return self._date

    # endregion

    ##
    # @brief Creates an instance of the class None date
    #
    # @return 
    def __init__(self):
        super().__init__()
        self._date = None

    ##
    # @brief Initialize the instance of the class with the given components.
    #
    # @param radLat Latitude, in radians. Negative value indicate latitude south.
    # @param radLon Longitude, in radians. Negative value indicate longitude west.
    # @param kmAlt Altitude above the ellipsoid model, in kilometers.
    # @param date The time associated with the coordinates.
    #
    # @return 
    def InitializeByRadLatAndRadLonAndKmAltAndDate(self, radLat, radLon, kmAlt, date):
        super().InitializeByRadLatAndRadLonAndKmAlt(radLat, radLon, kmAlt)
        self._date = date
        return self

    ##
    # @brief Initialize the instance of the class with the given Geo and Julian objects.
    #
    # @param geo The Geo object.
    # @param date The Julian date.
    #
    # @return 
    def InitializeByGeoAndDate(self, geo, date):
        super().InitializeByGeo(geo)
        self._date = date
        return self

    ##
    # @brief Initialize the instance of the class with the given ECI-time information.
    #
    # @param The ECI-time coordinates pair.
    #
    # @return 
    def InitializeByEci(self, eci):
        super().InitializeByEciAndDate(eci, eci.Date)
        self._date = eci.Date
        return self

    ##
    # @brief Initialize the instance of the class with the given ECI coordinates.
    #
    # @param eci The ECI coordinates.
    # @param date The Julian date.
    #
    # @return 
    def InitializeByEciAndDate(self, eci, date):
        super().InitializeByEciAndDate(eci, date)
        self._date = date
        return self



##
# @brief Class to encapsulate topo-centric coordinates.
class Topo(object):

    # region Properties

    @property
    ##
    # @brief The azimuth, in radians.
    #
    # @return 
    def AzimuthRad(self):
        return self._azimuthRad

    @property
    ##
    # @brief The elevation, in radians.
    #
    # @return 
    def ElevationRad(self):
        return self._elevationRad

    @property
    ##
    # @brief The azimuth, in degrees.
    #
    # @return 
    def AzimuthDeg(self):
        return Globals.ToDegrees(self._azimuthRad)

    @property
    ##
    # @brief The elevation, in degrees.
    #
    # @return 
    def ElevationDeg(self):
        return Globals.ToDegrees(self._elevationRad)

    @property
    ##
    # @brief The range, in kilometers.
    #
    # @return 
    def Range(self):
        return self._range

    @property
    ##
    # @brief The range rate, in kilometers per second.
    #        A negative value means "towards observer"
    #
    # @return 
    def RangeRate(self):
        return self._rangeRate

    # endregion

    ##
    # @brief Creates an instance of the class None AzimuthRad, None ElevationRad, None Range and 
    #        None RangeRate.
    #
    # @return 
    def __init__(self):
        self._azimuthRad    = None
        self._elevationRad  = None
        self._range         = None
        self._rangeRate     = None

    ##
    # @brief Initialize the instance of the class with the given componets.
    #
    # @param radAz Azimuth, in radians.
    # @param radEl Elevation, in radians.
    # @param r Range, in kilometers.
    # @param rangeRate Range rate, in kilometers per second. A negative
    #        range rate means "towards the observer".
    #
    # @return 
    def InitializeByRadAzAndRadElAndRangeAndRangeRate(self, radAz, radEl, r, rangeRate):
        self._azimuthRad    = radAz
        self._elevationRad  = radEl
        self._range         = r
        self._rangeRate     = rangeRate



##
# @brief Class to encapsulate topo-centric coordinates and a time.
class TopoTime(Topo):

    # region Properties

    @property
    ##
    # @brief The time associated with the coordinates.
    #
    # @return 
    def Date(self):
        return self._date

    # endregion

    ##
    # @brief Creates an instance of the class None date
    #
    # @return 
    def __init__(self):
        super().__init__()
        self._date = None

    ##
    # @brief Initialize the instance of the class with given topo and time information.
    #
    # @param topo
    # @param date
    #
    # @return 
    def InitializeByTopoAndDate(self, topo, date):
        super().InitializeByRadAzAndRadElAndRangeAndRangeRate(topo.AzimuthRad,
                                                              topo.ElevationRad,
                                                              topo.Range,
                                                              topo.RangeRate)
        this._date = date
        return self

    ##
    # @brief Initialize the instance of the class with given components.
    #
    # @param radAz Azimuth, in radians.
    # @param radEl Elevation, in radians.
    # @param r Range, in kilometers.
    # @param rangeRate Range rate, in kilometers per second. A negative
    #        range rate means "towards the observer".
    # @param date The time associated with the coordinates.
    #
    # @return 
    def InitializeByRadAzAndRadElAndRangeAndRangeRateAndDate(self,
                                                             radAz,
                                                             radEl,
                                                             r,
                                                             rangeRate,
                                                             date):
        super().InitializeByRadAzAndRadElAndRangeAndRangeRate(radAz,
                                                              radEl,
                                                              r,
                                                              rangeRate)
        self._date = date
        return self
