##
# @file NoradSGP4.py
# @brief 
# @author df_justforfun@163.com
# @version 1.0
# @date 2018-07-24

import math
from pythonOrbitTools.Orbit.NoradBase import NoradBase
# from pythonOrbitTools.Orbit.Orbit import Orbit
from pythonOrbitTools.Core.Globals import Globals

##
# @brief NORAD SGP4 implementation
class NoradSGP4(NoradBase):


    def __init__(self, orbit):
        super().__init__(orbit)

        self._m_c5      = 2.0 * self._m_coef1 * self._m_aodp * self._m_betao2 *\
                          (1.0 + 2.75 * (self._m_etasq + self._m_eeta) + self._m_eeta * self._m_etasq)
        self._m_omgcof  = self.Orbit.BStar * self._m_c3 * math.cos(self.Orbit.ArgPerigee)
        self._m_xmcof   = -(2.0 / 3.0) * self._m_coef * self.Orbit.BStar * Globals.Ae / self._m_eeta
        self._m_delmo   = math.pow(1.0 + self._m_eta * math.cos(self.Orbit.MeanAnomaly), 3.0)
        self._m_sinmo   = math.sin(self.Orbit.MeanAnomaly)

    ##
    # @brief Calculate satellite ECI position/velocity for a given time.
    #        This procedure returns the ECI position and velocity for the satellite
    #        in the orbit at the given number of minutes since the TLE epoch time.
    #        The algorithm uses NORAD`s Simplified General Perturbation 4 near earth
    #        orbit model.
    #
    # @param tsince Target time, in minutes-past-epoch format.
    #
    # @return AU-based position/velocity ECI coordinates.
    def GetPosition(self, tsince):
        # For m_perigee less than 220 kilometers, the isimp flag is set and
        # the equations are truncated to linear variation in square root of a
        # and quadratic variation in mean anomaly. Also, the m_c3 term, the
        # delta omega term, and the delta m term are dropped.

        isimp = False

        if (self._m_aodp * (1.0 - self._m_satEcc) / Globals.Ae) < (220.0 / Globals.Xkmper + Globals.Ae):
            isimp = True

        d2      = 0.0
        d3      = 0.0
        d4      = 0.0

        t3cof   = 0.0
        t4cof   = 0.0
        t5cof   = 0.0

        if not isimp:
            c1sq    = self._m_c1 * self._m_c1

            d2      = 4.0 * self._m_aodp * self._m_tsi * c1sq

            temp    = d2 * self._m_tsi * self._m_c1 / 3.0

            d3      = (17.0 * self._m_aodp + self._m_s4) * temp
            d4      = 0.5 * temp * self._m_aodp * self._m_tsi * (221.0 * self._m_aodp + 31.0 * self._m_s4) * self._m_c1
            t3cof   = d2 + 2.0 * c1sq
            t4cof   = 0.25 * (3.0 * d3 + self._m_c1 * (12.0 * d2 + 10.0 * c1sq))
            t5cof   = 0.2 * (3.0 * d4 + 12.0 * self._m_c1 * d3 + 6.0 * d2 * d2 + 15.0 * c1sq * (2.0 * d2 + c1sq))

        # Update for secular gravity and atmospheric drag.
        xmdf    = self.Orbit.MeanAnomaly + self._m_xmdot * tsince
        omgadf  = self.Orbit.ArgPerigee + self._m_omgdot * tsince
        xnoddf  = self.Orbit.RAAN + self._m_xnodot * tsince
        omega   = omgadf
        xmp     = xmdf
        tsq     = tsince * tsince
        xnode   = xnoddf + self._m_xnodcf * tsq
        tempa   = 1.0 - self._m_c1 * tsince
        tempe   = self.Orbit.BStar * self._m_c4 * tsince
        templ   = self._m_t2cof * tsq

        if not isimp:
            delomg  = self._m_omgcof * tsince
            delm    = self._m_xmcof * (math.pow(1.0 + self._m_eta * math.cos(xmdf), 3.0) - self._m_delmo)
            temp    = delomg + delm

            xmp     = xmdf + temp
            omega   = omgadf - temp

            tcube   = tsq * tsince
            tfour   = tsince * tcube

            tempa   = tempa - d2 * tsq - d3 * tcube - d4 * tfour
            tempe   = tempe + self.Orbit.BStar * self._m_c5 * (math.sin(xmp) - self._m_sinmo)
            templ   = templ + t3cof * tcube + tfour * (t4cof + tsince * t5cof)

        a   = self._m_aodp * Globals.Sqr(tempa)
        e   = self._m_satEcc - tempe

        xl  = xmp + omega + xnode + self._m_xnodp * templ
        xn  = Globals.Xke / math.pow(a, 1.5)

        return self.FinalPosition(self._m_satInc, omgadf, e, a, xl, xnode, xn, tsince)




