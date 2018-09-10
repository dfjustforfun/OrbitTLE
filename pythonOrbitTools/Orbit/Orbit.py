##
# @file Orbit.py
# @brief 
# @author df_justforfun@163.com
# @version 1.0
# @date 2018-07-18

import math
import datetime
from pythonOrbitTools.Core.Tle import Unit
from pythonOrbitTools.Core.Tle import Field
from pythonOrbitTools.Core.Globals import Globals
from pythonOrbitTools.Orbit.NoradSGP4 import NoradSGP4
from pythonOrbitTools.Orbit.NoradSDP4 import NoradSDP4

##
# @brief This class accepts a single satellite`s NORAD two-line element
#        set and provides information regarding the satellite`s orbit
#        such as period, axis length, ECI coordinates, velocity, etc.
class Orbit(object):

    # region Properties

    @property
    def Tle(self):
        return self._tle

    @property
    def TleLine1(self):
        return self._tle.Line1

    @property
    def TleLine2(self):
        return self._tle.Line2

    @property
    def Epoch(self):
        return self._epoch # Julian

    @property
    def EpochTime(self):
        return self._epoch.ToTime() # datetime

    @property
    def NoradModel(self):
        return self._noradModel

    # region "Recovered" from the input elements

    @property
    def SemiMajor(self):
        return self._m_aeAxisSemiMajorRec

    @property
    def SemiMinor(self):
        return self._m_aeAxisSemiMinorRec

    @property
    def MeanMotion(self):
        return self._m_rmMeanMotionRec

    @property
    def Major(self):
        return 2.0*self.SemiMajor

    @property
    def Minor(self):
        return 2.0*self.SemiMinor

    @property
    def Perigee(self):
        return self._m_kmPerigeeRec

    @property
    def Apogee(self):
        return self._m_kmApogeeRec

    # endregion

    @property
    def Inclination(self):
        return self._m_Inclination

    @property
    def Eccentricity(self):
        return self._m_Eccentricity

    @property
    def RAAN(self):
        return self._m_RAAN

    @property
    def ArgPerigee(self):
        return self._m_ArgPerigee

    @property
    def BStar(self):
        return self._m_BStar

    @property
    def Drag(self):
        return self._m_Drag

    @property
    def MeanAnomaly(self):
        return self._m_MeanAnomaly

    @property
    def TleMeanMotion(self):
        return self._m_TleMeanMotion

    @property
    def SatNoradId(self):
        return self._tle.NoradNumber

    @property
    def SatName(self):
        return self._tle.Name

    @property
    def SatNameLong(self):
        return self.SatName + " #" + self.SatNoradId

    @property
    def Period(self):
        if self._m_Period.total_seconds() < 0.0:
            # Calculate the period using the recovered mean motion.
            if self.MeanMotion == 0.0:
                self._m_Period = datetime.timedelta()
            else:
                secs = (Globals.TwoPi/self.MeanMotion)*60.0
                self._m_Period = datetime.timedelta(seconds=secs)

        return self._m_Period

    # endregion


    ##
    # @brief Standard constructor
    #
    # @param tle Two-line element orbital parameters.
    #
    # @return 
    def __init__(self, tle):
        self._tle   = tle
        self._epoch = tle.EpochJulian

        # Caching variables
        self._m_Period = datetime.timedelta(seconds=-1)

        # TLE caching variables
        self._m_Inclination     = self.GetRad(Field.Inclination)
        self._m_Eccentricity    = self._tle.GetFieldAsValue(Field.Eccentricity)
        self._m_RAAN            = self.GetRad(Field.Raan)
        self._m_ArgPerigee      = self.GetRad(Field.ArgPerigee)
        self._m_BStar           = self._tle.GetFieldAsValue(Field.BStarDrag)
        self._m_Drag            = self._tle.GetFieldAsValue(Field.MeanMotionDt)
        self._m_MeanAnomaly     = self.GetRad(Field.MeanAnomaly)
        self._m_TleMeanMotion   = self._tle.GetFieldAsValue(Field.MeanMotion)

        # Recover the original mean motion and semimajor axis from the input elements.
        mm      = self.TleMeanMotion
        rpmin   = mm*Globals.TwoPi/Globals.MinPerDay # rads per minute

        a1      = math.pow(Globals.Xke/rpmin, 2.0/3.0)
        e       = self.Eccentricity
        i       = self.Inclination
        temp    = (1.5*Globals.Ck2*(3.0*Globals.Sqr(math.cos(i))-1.0)/math.pow(1.0-e*e, 1.5))
        delta1  = temp/(a1*a1)
        a0      = a1*(1.0-delta1*((1.0/3.0)+delta1*(1.0+134.0/81.0*delta1)))

        delta0  = temp/(a0*a0)

        # Caching variables recovered from the input TLE elements
        self._m_rmMeanMotionRec     = rpmin/(1.0+delta0) # radians per minute
        self._m_aeAxisSemiMajorRec  = a0/(1.0-delta0) # semimajor axis, in AE units
        self._m_aeAxisSemiMinorRec  = self._m_aeAxisSemiMajorRec*math.sqrt(1.0-(e*e)) # semiminor axis, in AE units.
        self._m_kmPerigeeRec        = Globals.Xkmper*(self._m_aeAxisSemiMajorRec*(1.0-e)-Globals.Ae) # perigee, in km
        self._m_kmApogeeRec         = Globals.Xkmper*(self._m_aeAxisSemiMajorRec*(1.0+e)-Globals.Ae) # apogee, in km


        if (self.Period.total_seconds() / 60.0) >= 225.0:
            # SDP4 - period >= 225 minutes.
            self._noradModel = NoradSDP4(self)
        else:
            # SGP4 - period < 225 minutes.
            self._noradModel = NoradSGP4(self)

    ##
    # @brief Calculate ECI position/velocity for a given time.
    #
    # @param mpe Target time, in minutes past the TLE epoch.
    #
    # @return Kilometer-based position/velocity ECI coordinates. 
    def PositionEciByMpe(self, mpe):

        eci         = self.NoradModel.GetPosition(mpe)

        # Convert ECI vector units from AU to kilometers.
        radiusAe    = Globals.Xkmper / Globals.Ae

        eci.ScalePosVector(radiusAe) # km
        eci.ScaleVelVector(radiusAe * (Globals.MinPerDay / 86400.0)) # km /sec

        return eci

    ##
    # @brief Calculate ECI position/velocity for a given time.
    #
    # @param utc Target time (UTC).
    #
    # @return Kilometer-based position/velocity ECI coordinates. 
    def PositionEciByDateTime(self, utc):
        return self.PositionEciByMpe(self.TPlusEpoch(utc).total_seconds()/60.0)


    ##
    # @brief Returns elspaed time from epoch to given time or current time. 
    #        "Predicted" TLEs can have epochs in the future.
    #
    # @param utc 
    #
    # @return 
    def TPlusEpoch(self, utc=datetime.datetime.utcnow()):
        return utc - self.EpochTime


    # region Utility

    def GetRad(self, fld):
        return self._tle.GetFieldAsValue(fld, Unit.Radians)

    def GetDeg(self, fld):
        return self._tle.GetFieldAsValue(fld, Unit.Degrees)

    # endregion
