##
# @file NoradBase.py
# @brief 
# @author df_justforfun@163.com
# @version 1.0
# @date 2018-07-23

import math
import datetime
from abc import ABCMeta, abstractmethod
# from pythonOrbitTools.Orbit.Orbit import Orbit
from pythonOrbitTools.Core.Globals import Globals
from pythonOrbitTools.Core.Vector import Vector
from pythonOrbitTools.Core.Eci import EciTime
from pythonOrbitTools.Core.Julian import Julian

##
# @brief This class provides a base class for the NORAD SGP4/SDP4 orbit models.
class NoradBase(metaclass=ABCMeta):

    @property
    def Orbit(self):
        return self._orbit

    @abstractmethod
    def GetPosition(self, tsince):
        pass

    def __init__(self, orbit):
        self._orbit = orbit
        self.Initialize()

    ##
    # @brief Perform the initialization of member variables, specifically the variables
    #        used by derived-class objects to calculate ECI coordinates.
    #
    # @return 
    def Initialize(self):
        # Orbital parameter variables which need only the calculated one
        # time for a given orbit (ECI position time-independent).
        # Initialize any variables which are time-independent when
        # calculating the ECI coordinates of the satellite.

        self._m_satInc      = self.Orbit.Inclination # inclination
        self._m_satEcc      = self.Orbit.Eccentricity # eccentricity

        self._m_cosio       = math.cos(self._m_satInc)
        self._m_theta2      = self._m_cosio * self._m_cosio
        self._m_x3thm1      = 3.0 * self._m_theta2 - 1.0
        self._m_eosq        = self._m_satEcc * self._m_satEcc
        self._m_betao2      = 1.0 - self._m_eosq
        self._m_betao       = math.sqrt(self._m_betao2)

        # The "recovered" semimajor axis and mean motion.
        self._m_aodp        = self.Orbit.SemiMajor
        self._m_xnodp       = self.Orbit.MeanMotion

        # For perigee below 156 km, the values of Globals.S and Globals.QQMS2T are altered.
        self._m_perigee     = Globals.Xkmper * (self._m_aodp * (1.0 - self._m_satEcc) - Globals.Ae)

        self._m_s4          = Globals.S
        self._m_qoms24      = Globals.Qoms2t

        if self._m_perigee < 156.0:
            self._m_s4 = self._m_perigee - 78.0

            if self._m_perigee <= 98.0:
                self._m_s4 = 20.0

            self._m_qoms24 = math.pow((120.0 - self._m_s4) * Globals.Ae / Globals.Xkmper, 4.0)
            self._m_s4     = self._m_s4 / Globals.Xkmper + Globals.Ae

        pinvsq = 1.0 / (self._m_aodp * self._m_aodp * self._m_betao2 * self._m_betao2)

        self._m_tsi     = 1.0 / (self._m_aodp - self._m_s4)
        self._m_eta     = self._m_aodp * self._m_satEcc * self._m_tsi
        self._m_etasq   = self._m_eta * self._m_eta
        self._m_eeta    = self._m_satEcc * self._m_eta

        psisq = math.fabs(1.0 - self._m_etasq)

        self._m_coef    = self._m_qoms24 * math.pow(self._m_tsi, 4.0)
        self._m_coef1   = self._m_coef / math.pow(psisq, 3.5)

        c2 = self._m_coef1 * self._m_xnodp *\
                (self._m_aodp * (1.0 + 1.5 * self._m_etasq + self._m_eeta * (4.0 + self._m_etasq)) +\
                 0.75 * Globals.Ck2 * self._m_tsi / psisq * self._m_x3thm1 *\
                 (0.8 + 3.0 * self._m_etasq * (8.0 + self._m_etasq)))


        self._m_c1      = self.Orbit.BStar * c2
        self._m_sinio   = math.sin(self._m_satInc)

        a3ovk2 = -Globals.Xj3 / Globals.Ck2 * math.pow(Globals.Ae, 3.0)

        self._m_c3      = self._m_coef * self._m_tsi * a3ovk2 * self._m_xnodp * Globals.Ae * self._m_sinio / self._m_satEcc
        self._m_x1mth2  = 1.0 - self._m_theta2
        self._m_c4      = 2.0 * self._m_xnodp * self._m_coef1 * self._m_aodp * self._m_betao2 *\
                          (self._m_eta * (2.0 + 0.5 * self._m_etasq) +\
                           self._m_satEcc * (0.5 + 2.0 * self._m_etasq) -\
                           2.0 * Globals.Ck2 * self._m_tsi / (self._m_aodp * psisq) *\
                           (-3.0 * self._m_x3thm1 * (1.0 - 2.0 * self._m_eeta + self._m_etasq * (1.5 -0.5 * self._m_eeta)) +
                            0.75 * self._m_x1mth2 *\
                            (2.0 * self._m_etasq - self._m_eeta * (1.0 + self._m_etasq)) *\
                            math.cos(2.0 * self.Orbit.ArgPerigee)))

        theta4 = self._m_theta2 * self._m_theta2
        temp1  = 3.0 * Globals.Ck2 * pinvsq * self._m_xnodp
        temp2  = temp1 * Globals.Ck2 * pinvsq
        temp3  = 1.25 * Globals.Ck4 * pinvsq * pinvsq * self._m_xnodp

        self._m_xmdot   = self._m_xnodp + 0.5 * temp1 * self._m_betao * self._m_x3thm1 +\
                          0.0625 * temp2 * self._m_betao *\
                          (13.0 - 78.0 * self._m_theta2 + 137.0 * theta4)

        x1m5th = 1.0 - 5.0 * self._m_theta2

        self._m_omgdot  = -0.5 * temp1 * x1m5th + 0.0625 * temp2 *\
                          (7.0 - 114.0 * self._m_theta2 + 395.0 * theta4) +\
                          temp3 * (3.0 - 36.0 * self._m_theta2 + 49.0 * theta4)

        xhdot1 = -temp1 * self._m_cosio

        self._m_xnodot  = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * self._m_theta2) +\
                          2.0 * temp3 * (3.0 - 7.0 * self._m_theta2)) * self._m_cosio
        self._m_xnodcf  = 3.5 * self._m_betao2 * xhdot1 * self._m_c1
        self._m_t2cof   = 1.5 * self._m_c1
        self._m_xlcof   = 0.125 * a3ovk2 * self._m_sinio * (3.0 + 5.0 * self._m_cosio) / (1.0 + self._m_cosio)
        self._m_aycof   = 0.25 * a3ovk2 * self._m_sinio
        self._m_x7thm1  = 7.0 * self._m_theta2 - 1.0


    def FinalPosition(self, incl, omega, e, a, xl, xnode, xn, tsince):
        if (e*e) > 1.0:
            raise ValueError("Error in satellite data")

        beta = math.sqrt(1.0 - e*e)

        # Long period periodics 
        axn     = e * math.cos(omega)
        temp    = 1.0 / (a * beta * beta)
        xll     = temp * self._m_xlcof * axn
        aynl    = temp * self._m_aycof
        xlt     = xl + xll
        ayn     = e * math.sin(omega) + aynl

        # Solve Kepler`s Equation
        capu    = Globals.Fmod2p(xlt - xnode)
        temp2   = capu
        temp3   = 0.0
        temp4   = 0.0
        temp5   = 0.0
        temp6   = 0.0
        sinepw  = 0.0
        cosepw  = 0.0
        fDone   = False

        i = 1
        while ((i <= 10) and (not fDone)):
            sinepw  = math.sin(temp2)
            cosepw  = math.cos(temp2)
            temp3   = axn * sinepw
            temp4   = ayn * cosepw
            temp5   = axn * cosepw
            temp6   = ayn * sinepw

            epw     = (capu - temp4 + temp3 - temp2) / (1.0 - temp5 -temp6) + temp2

            if math.fabs(epw - temp2) <= 1.0e-06:
                fDone = True
            else:
                temp2 = epw

            i += 1

        # Short period preliminary quantities
        ecose   = temp5 + temp6
        esine   = temp3 - temp4
        elsq    = axn * axn + ayn * ayn
        temp    = 1.0 - elsq
        pl      = a * temp
        r       = a * (1.0 - ecose)
        temp1   = 1.0 / r
        rdot    = Globals.Xke * math.sqrt(a) * esine * temp1
        rfdot   = Globals.Xke * math.sqrt(pl) * temp1
        temp2   = a * temp1
        betal   = math.sqrt(temp)
        temp3   = 1.0 / (1.0 + betal)
        cosu    = temp2 * (cosepw - axn + ayn * esine * temp3)
        sinu    = temp2 * (sinepw - ayn - axn * esine * temp3)
        u       = Globals.AcTan(sinu, cosu)
        sin2u   = 2.0 * sinu * cosu
        cos2u   = 2.0 * cosu * cosu - 1.0


        temp    = 1.0 / pl
        temp1   = Globals.Ck2 * temp
        temp2   = temp1 * temp

        # Update for short periodics
        rk      = r * (1.0 - 1.5 * temp2 * betal * self._m_x3thm1) + 0.5 * temp1 * self._m_x1mth2 * cos2u
        uk      = u - 0.25 * temp2 * self._m_x7thm1 * sin2u
        xnodek  = xnode + 1.5 * temp2 * self._m_cosio * sin2u
        xinck   = incl + 1.5 * temp2 * self._m_cosio * self._m_sinio * cos2u
        rdotk   = rdot - xn * temp1 * self._m_x1mth2 * sin2u
        rfdotk  = rfdot + xn * temp1 * (self._m_x1mth2 * cos2u + 1.5 * self._m_x3thm1)

        # Orientation vectors
        sinuk   = math.sin(uk)
        cosuk   = math.cos(uk)
        sinik   = math.sin(xinck)
        cosik   = math.cos(xinck)
        sinnok  = math.sin(xnodek)
        cosnok  = math.cos(xnodek)
        xmx     = -sinnok * cosik
        xmy     =  cosnok * cosik
        ux      = xmx * sinuk + cosnok * cosuk
        uy      = xmy * sinuk + sinnok * cosuk
        uz      = sinik * sinuk
        vx      = xmx * cosuk - cosnok * sinuk
        vy      = xmy * cosuk - sinnok * sinuk
        vz      = sinik * cosuk

        # Position
        x       = rk * ux
        y       = rk * uy
        z       = rk * uz

        vecPos  = Vector(x, y, z)
        gmt     = self.Orbit.EpochTime + datetime.timedelta(minutes=tsince)

        # Validate on altitude
        altKm   = (vecPos.Magnitude() * (Globals.Xkmper / Globals.Ae))
        if altKm < Globals.Xkmper:
            raise ValueError(str(gmt)+self.Orbit.SatNameLong)

        # Velocity
        xdot    = rdotk * ux + rfdotk * vx
        ydot    = rdotk * uy + rfdotk * vy
        zdot    = rdotk * uz + rfdotk * vz

        vecVel  = Vector(xdot, ydot, zdot)

        return EciTime().InitializeByPosAndVelAndDate(vecPos, vecVel, Julian().InitializeByUTC(gmt))


if __name__ == "__main__":
    pass

