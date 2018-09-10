##
# @file NoradSDP4.py
# @brief 
# @author df_justforfun@163.com
# @version 1.0
# @date 2018-08-17

import math
from pythonOrbitTools.Orbit.NoradBase import NoradBase
from pythonOrbitTools.Core.Globals import Globals

##
# @brief NORAD SDP4 implementation
class NoradSDP4(NoradBase):

    # region Constants

    zns     = 1.19459E-5
    zes     = 0.01675
    znl     = 1.5835218E-4
    zel     = 0.05490
    thdt    = 4.3752691E-3

    # endregion

    def __init__(self, orbit):
        super().__init__(orbit)

        sinarg = math.sin(self.Orbit.ArgPerigee)
        cosarg = math.cos(self.Orbit.ArgPerigee)

        # Deep space initialization
        jd = self.Orbit.Epoch

        self._dp_thgr = jd.ToGmst()

        eq      = self.Orbit.Eccentricity
        aqnv    = 1.0 / self.Orbit.SemiMajor

        self._dp_xqncl = self.Orbit.Inclination

        xmao    = self.Orbit.MeanAnomaly
        xpidot  = self._m_omgdot + self._m_xnodot
        sinq    = math.sin(self.Orbit.RAAN)
        cosq    = math.cos(self.Orbit.RAAN)

        self._dp_omegaq = self.Orbit.ArgPerigee

        # region Lunar / Solar terms

        # Initialize lunar solar terms
        day = jd.FromJan0_12h_1900

        dpi_xnodce    = 4.5236020 - 9.2422029E-4 * day
        dpi_stem      = math.sin(dpi_xnodce)
        dpi_ctem      = math.cos(dpi_xnodce)
        dpi_zcosil    = 0.91375164 - 0.03568096 * dpi_ctem
        dpi_zsinil    = math.sqrt(1.0 - dpi_zcosil * dpi_zcosil)
        dpi_zsinhl    = 0.089683511 * dpi_stem / dpi_zsinil
        dpi_zcoshl    = math.sqrt(1.0 - dpi_zsinhl * dpi_zsinhl)
        dpi_c         = 4.7199672 + 0.22997150 * day
        dpi_gam       = 5.8351514 + 0.0019443680 * day

        self._dp_zmol = Globals.Fmod2p(dpi_c - dpi_gam)

        dpi_zx = 0.39785416 * dpi_stem / dpi_zsinil
        dpi_zy = dpi_zcoshl * dpi_ctem + 0.91744867 * dpi_zsinhl * dpi_stem

        dpi_zx = Globals.AcTan(dpi_zx, dpi_zy) + dpi_gam - dpi_xnodce

        dpi_zcosgl = math.cos(dpi_zx)
        dpi_zsingl = math.sin(dpi_zx)

        dp_zmos = 6.2565837 + 0.017201977 * day
        dp_zmos = Globals.Fmod2p(dp_zmos)

        zcosis  =  0.91744867
        zsinis  =  0.39785416
        zsings  = -0.98088458
        zcosgs  =  0.1945905
        c1ss    =  2.9864797E-6

        zcosg   = zcosgs
        zsing   = zsings
        zcosi   = zcosis
        zsini   = zsinis
        zcosh   = cosq
        zsinh   = sinq
        cc      = c1ss
        zn      = NoradSDP4.zns
        ze      = NoradSDP4.zes
        xnoi    = 1.0 / self.Orbit.MeanMotion

        se      = 0.0
        si      = 0.0
        sl      = 0.0
        sgh     = 0.0
        sh      = 0.0

        eosq    = Globals.Sqr(self.Orbit.Eccentricity)

        # Apply the solar and lunar terms on the first pass, then re-apply the
        # solar terms again on the second pass

        for i in range(1, 3):
            # Do solar terms
            a1  =  zcosg * zcosh + zsing * zcosi * zsinh
            a3  = -zsing * zcosh + zcosg * zcosi * zsinh
            a7  = -zcosg * zsinh + zsing * zcosi * zcosh
            a8  =  zsing * zsini
            a9  =  zsing * zsinh + zcosg * zcosi * zcosh
            a10 =  zcosg * zsini
            a2  =  self._m_cosio * a7 + self._m_sinio * a8
            a4  =  self._m_cosio * a9 + self._m_sinio * a10
            a5  = -self._m_sinio * a7 + self._m_cosio * a8
            a6  = -self._m_sinio * a9 + self._m_cosio * a10
            x1  =  a1 * cosarg + a2 * sinarg
            x2  =  a3 * cosarg + a4 * sinarg
            x3  = -a1 * sinarg + a2 * cosarg
            x4  = -a3 * sinarg + a4 * cosarg
            x5  =  a5 * sinarg
            x6  =  a6 * sinarg
            x7  =  a5 * cosarg
            x8  =  a6 * cosarg
            z31 =  12.0 * x1 * x1 - 3.0 * x3 * x3
            z32 =  24.0 * x1 * x1 - 6.0 * x3 * x4
            z33 =  12.0 * x2 * x2 - 3.0 * x4 * x4
            z1  =  3.0 * (a1 * a1 + a2 * a2) + z31 * eosq
            z2  =  6.0 * (a1 * a3 + a2 * a4) + z32 * eosq
            z3  =  3.0 * (a3 * a3 + a4 * a4) + z33 * eosq
            z11 = -6.0 * a1 * a5 + eosq * (-24.0 * x1 * x7 - 6.0 * x3 * x5)
            z12 = -6.0 * (a1 * a6 + a3 * a5) + \
                    eosq * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5))
            z13 = -6.0 * a3 * a6 + eosq * (-24.0 * x2 * x8 - 6.0 * x4 * x6)
            z21 =  6.0 * a2 * a5 + eosq * (24.0 * x1 * x5 - 6.0 * x3 * x7)
            z22 =  6.0 * (a4 * a5 + a2 * a6) + \
                    eosq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8))
            z23 =  6.0 * a4 * a6 + eosq * (24.0 * x2 * x6 - 6.0 * x4 * x8)
            z1  =  z1 + z1 + self._m_betao2 * z31
            z2  =  z2 + z2 + self._m_betao2 * z32
            z3  =  z3 + z3 + self._m_betao2 * z33
            s3  =  cc * xnoi
            s2  = -0.5 * s3 / self._m_betao
            s4  =  s3 * self._m_betao
            s1  = -15.0 * eq * s4
            s5  =  x1 * x3 + x2 * x4
            s6  =  x2 * x3 + x1 * x4
            s7  =  x2 * x4 - x1 * x3
            se  =  s1 * zn * s5
            si  =  s2 * zn * (z11 + z13)
            sl  = -zn * s3 * (z1 + z3 - 14.0 - 6.0 * eosq)
            sgh =  s4 * zn * (z31 + z33 - 6.0)

            if self.Orbit.Inclination < 5.2359877E-2:
                sh = 0.0
            else:
                sh = -zn * s2 * (z21 + z23)

            self._dp_ee2    =  2.0 * s1 * s6
            self._dp_e3     =  2.0 * s1 * s7
            self._dp_xi2    =  2.0 * s2 * z12
            self._dp_xi3    =  2.0 * s2 * (z13 - z11)
            self._dp_xl2    = -2.0 * s3 * z2
            self._dp_xl3    = -2.0 * s3 * (z3 - z1)
            self._dp_xl4    = -2.0 * s3 * (-21.0 - 9.0 * eosq) * ze
            self._dp_xgh2   =  2.0 * s4 * z32
            self._dp_xgh3   =  2.0 * s4 * (z33 - z31)
            self._dp_xgh4   = -18.0 * s4 * ze
            self._dp_xh2    = -2.0 * s2 * z22
            self._dp_xh3    = -2.0 * s2 * (z23 - z21)

            if i == 1:
                # Do lunar terms
                self._dp_sse    = se
                self._dp_ssi    = si
                self._dp_ssl    = sl
                self._dp_ssh    = sh / self._m_sinio
                self._dp_ssg    = sgh - self._m_cosio * self._dp_ssh
                self._dp_se2    = self._dp_ee2
                self._dp_si2    = self._dp_xi2
                self._dp_sl2    = self._dp_xl2
                self._dp_sgh2   = self._dp_xgh2
                self._dp_sh2    = self._dp_xh2
                self._dp_se3    = self._dp_e3
                self._dp_si3    = self._dp_xi3
                self._dp_sl3    = self._dp_xl3
                self._dp_sgh3   = self._dp_xgh3
                self._dp_sh3    = self._dp_xh3
                self._dp_sl4    = self._dp_xl4
                self._dp_sgh4   = self._dp_xgh4
                zcosg           = dpi_zcosgl
                zsing           = dpi_zsingl
                zcosi           = dpi_zcosil
                zsini           = dpi_zsinil
                zcosh           = dpi_zcoshl * cosq + dpi_zsinhl * sinq
                zsinh           = sinq * dpi_zcoshl - cosq * dpi_zsinhl
                zn              = NoradSDP4.znl

                c1l             = 4.7968065E-7

                cc              = c1l
                ze              = NoradSDP4.zel


        # endregion
        self._dp_sse = self._dp_sse + se
        self._dp_ssi = self._dp_ssi + si
        self._dp_ssl = self._dp_ssl + sl
        self._dp_ssg = self._dp_ssg + sgh - self._m_cosio / self._m_sinio * sh
        self._dp_ssh = self._dp_ssh + sh / self._m_sinio

        # Geopotential resonance initialization for 12 hour orbits
        self._gp_reso = False # geopotential resonant
        self._gp_sync = False # geopotential synchronous

        g310    = 0.0
        f220    = 0.0
        bfact   = 0.0

        # Determine if orbit is 12- or 24-hour resonant.
        # Mean motion is given in radians per minute.
        if self.Orbit.MeanMotion > 0.0034906585 and self.Orbit.MeanMotion < 0.0052359877:
            # Orbit is within the Clarke Belt (period is 24-hour resonant)
            # Synchronous reasonance terms initialization
            self._gp_reso = True
            self._gp_sync = True

            # region 24-hour resonant

            g200    = 1.0 + eosq * (-2.5 + 0.8125 * eosq)
            g310    = 1.0 + 2.0 * eosq
            g300    = 1.0 + eosq * (-6.0 + 6.60937 * eosq)
            f220    = 0.75 * (1.0 + self._m_cosio) * (1.0 + self._m_cosio)
            f311    = 0.9375 * self._m_sinio * self._m_sinio * (1.0 + 3 * self._m_cosio) - 0.75 * (1.0 + self._m_cosio)
            f330    = 1.0 + self._m_cosio
            f330    = 1.875 * f330 * f330 * f330

            q22     = 1.7891679e-06
            q33     = 2.2123015e-07
            q31     = 2.1460748e-06

            self._dp_del1   = 3.0 * self._m_xnodp * self._m_xnodp * aqnv * aqnv
            self._dp_del2   = 2.0 * self._dp_del1 * f220 * g200 * q22
            self._dp_del3   = 3.0 * self._dp_del1 * f330 * g300 * q33 * aqnv
            self._dp_del1   = self._dp_del1 * f311 * g310 * q31 * aqnv
            self._dp_xlamo  = xmao + self.Orbit.RAAN + self.Orbit.ArgPerigee - self._dp_thgr
            bfact           = self._m_xmdot + xpidot - NoradSDP4.thdt
            bfact           = bfact + self._dp_ssl + self._dp_ssg + self._dp_ssh

            # endregion
        elif self.Orbit.MeanMotion >= 8.26E-3 and self.Orbit.MeanMotion <= 9.24E-3 and eq >= 0.5:
            # Period is 12-hour resonant
            self._gp_reso = True

            # region 12-hour resonant

            eoc     =  eq * eosq
            g201    = -0.306 - (eq - 0.64) * 0.440

            if eq <= 0.65:
                g211 =  3.616 - 13.247 * eq + 16.290 * eosq
                g310 = -19.302 + 117.390 * eq - 228.419 * eosq + 156.591 * eoc
                g322 = -18.9068 + 109.7927 * eq - 214.6334 * eosq + 146.5816 * eoc
                g410 = -41.122 + 242.694 * eq - 471.094 * eosq + 313.953 * eoc
                g422 = -146.407 + 841.880 * eq - 1629.014 * eosq + 1083.435 * eoc
                g520 = -532.114 + 3017.977 * eq - 5740.0 * eosq + 3708.276 * eoc
            else:
                g211 = -72.099 + 331.819 * eq - 508.738 * eosq + 266.724 * eoc
                g310 = -346.844 + 1582.851 * eq - 2415.925 * eosq + 1246.113 * eoc
                g322 = -342.585 + 1554.908 * eq - 2366.899 * eosq + 1215.972 * eoc
                g410 = -1052.797 + 4758.686 * eq - 7193.992 * eosq + 3651.957 * eoc
                g422 = -3581.69 + 16178.11 * eq - 24462.77 * eosq + 12422.52 * eoc

                if eq <= 0.715:
                    g520 = 1464.74 - 4664.75 * eq + 3763.64 * eosq
                else:
                    g520 = -5149.66 + 29936.92 * eq - 54087.36 * eosq + 31324.56 * eoc

            if eq < 0.7:
                g533 = -919.2277 + 4988.61 * eq - 9064.77 * eosq + 5542.21 * eoc
                g521 = -822.71072 + 4568.6173 * eq - 8491.4146 * eosq + 5337.524 * eoc
                g532 = -853.666 + 4690.25 * eq - 8624.77 * eosq + 5341.4 * eoc
            else:
                g533 = -37995.78 + 161616.52 * eq - 229838.2 * eosq + 109377.94 * eoc
                g521 = -51752.104 + 218913.95 * eq - 309468.16 * eosq + 146349.42 * eoc
                g532 = -40023.88 + 170470.89 * eq - 242699.48 * eosq + 115605.82 * eoc

            sini2 = self._m_sinio * self._m_sinio
            cosi2 = self._m_cosio * self._m_cosio

            f220    =  0.75 * (1.0 + 2.0 * self._m_cosio + cosi2)
            f221    =  1.5 * sini2
            f321    =  1.875 * self._m_sinio * (1.0 - 2.0 * self._m_cosio - 3.0 * cosi2)
            f322    = -1.875 * self._m_sinio * (1.0 + 2.0 * self._m_cosio - 3.0 * cosi2)
            f441    =  35.0 * sini2 * f220
            f442    =  39.3750 * sini2 * sini2
            f522    =  9.84375 * self._m_sinio * (sini2 * (1.0 - 2.0 * self._m_cosio - 5.0 * cosi2) + 0.33333333 * (-2.0 + 4.0 * self._m_cosio + 6.0 * cosi2))
            f523    =  self._m_sinio * (4.92187512 * sini2 * (-2.0 - 4.0 * self._m_cosio + 10.0 * cosi2) + 6.56250012 * (1.0 + 2.0 * self._m_cosio - 3.0 * cosi2))
            f542    =  29.53125 * self._m_sinio * ( 2.0 - 8.0 * self._m_cosio + cosi2 * (-12.0 + 8.0 * self._m_cosio + 10.0 * cosi2))
            f543    =  29.53125 * self._m_sinio * (-2.0 - 8.0 * self._m_cosio + cosi2 * ( 12.0 + 8.0 * self._m_cosio - 10.0 * cosi2))
            xno2    =  self._m_xnodp * self._m_xnodp
            ainv2   =  aqnv * aqnv
            temp1   =  3.0 * xno2 * ainv2

            root22  = 1.7891679E-6
            root32  = 3.7393792E-7
            root44  = 7.3636953E-9
            root52  = 1.1428639E-7
            root54  = 2.1765803E-9

            temp    = temp1 * root22

            self._dp_d2201  = temp * f220 * g201
            self._dp_d2211  = temp * f221 * g211
            temp1           = temp1 * aqnv
            temp            = temp1 * root32
            self._dp_d3210  = temp * f321 * g310
            self._dp_d3222  = temp * f322 * g322
            temp1           = temp1 * aqnv
            temp            = 2.0 * temp1 * root44
            self._dp_d4410  = temp * f441 * g410
            self._dp_d4422  = temp * f442 * g422
            temp1           = temp1 * aqnv
            temp            = temp1 * root52
            self._dp_d5220  = temp * f522 * g520
            self._dp_d5232  = temp * f523 * g532
            temp            = 2.0 * temp1 * root54
            self._dp_d5421  = temp * f542 * g521
            self._dp_d5433  = temp * f543 * g533
            self._dp_xlamo  = xmao + self.Orbit.RAAN + self.Orbit.RAAN - self._dp_thgr - self._dp_thgr
            bfact           = self._m_xmdot + self._m_xnodot + self._m_xnodot - NoradSDP4.thdt -NoradSDP4.thdt
            bfact           = bfact + self._dp_ssl + self._dp_ssh + self._dp_ssh

            # endregion
        if self._gp_reso or self._gp_sync:
            self._dp_xfact  = bfact - self._m_xnodp

            # Initialize integrator
            self._dp_xli    =  self._dp_xlamo
            self._dp_xni    =  self._m_xnodp
            dp_atime        =  0.0; # performed by runtime
            self._dp_stepp  =  720.0
            self._dp_stepn  = -720.0
            self._dp_step2  =  259200.0

    def DeepCalcDotTerms(self, pxndot, pxnddt, pxldot):
        # Dot terms calculated
        if self._gp_sync:
            fasx2   = 0.13130908
            fasx4   = 2.8843198
            fasx6   = 0.37448087

            pxndot  = self._dp_del1 * math.sin(self._dp_xli - fasx2) + \
                        self._dp_del2 * math.sin(2.0 * (self._dp_xli - fasx4)) + \
                        self._dp_del3 * math.sin(3.0 * (self._dp_xli - fasx6))
            pxnddt  = self._dp_del1 * math.cos(self._dp_xli - fasx2) + \
                        2.0 * self._dp_del2 * math.cos(2.0 * (self._dp_xli - fasx4)) + \
                        3.0 * self._dp_del3 * math.cos(3.0 * (self._dp_xli - fasx6))
        else:
            g54     = 4.4108898
            g52     = 1.0508330
            g44     = 1.8014998
            g22     = 5.7686396
            g32     = 0.95240898

            xomi    = self._dp_omegaq + self._m_omgdot * self._dp_atime
            x2omi   = xomi + xomi
            x2li    = self._dp_xli + self._dp_xli

            pxndot  = self._dp_d2201 * math.sin(x2omi + self._dp_xli - g22) + \
                        self._dp_d2211 * math.sin(self._dp_xli - g22) + \
                        self._dp_d3210 * math.sin( xomi + self._dp_xli - g32) + \
                        self._dp_d3222 * math.sin(-xomi + self._dp_xli - g32) + \
                        self._dp_d4410 * math.sin(x2omi + x2li - g44) + \
                        self._dp_d4422 * math.sin(x2li - g44) + \
                        self._dp_d5220 * math.sin( xomi + self._dp_xli - g52) + \
                        self._dp_d5232 * math.sin(-xomi + self._dp_xli - g52) + \
                        self._dp_d5421 * math.sin( xomi + x2li - g54) + \
                        self._dp_d5433 * math.sin(-xomi + x2li - g54)

            pxnddt  = self._dp_d2201 * math.cos(x2omi + self._dp_xli - g22) + \
                        self._dp_d2211 * math.cos(self._dp_xli - g22) + \
                        self._dp_d3210 * math.cos( xomi + self._dp_xli - g32) + \
                        self._dp_d3222 * math.cos(-xomi + self._dp_xli - g32) + \
                        self._dp_d5220 * math.cos( xomi + self._dp_xli - g52) + \
                        self._dp_d5232 * math.cos(-xomi + self._dp_xli - g52) + \
                        2.0 * (self._dp_d4410 * math.cos(x2omi + x2li - g44) + \
                        self._dp_d4422 * math.cos(x2li - g44) + \
                        self._dp_d5421 * math.cos( xomi + x2li - g54) + \
                        self._dp_d5433 * math.cos(-xomi + x2li - g54))

        pxldot = self._dp_xni + self._dp_xfact
        pxnddt = pxnddt * pxldot

        return True, pxndot, pxnddt, pxldot

    def DeepCalcIntegrator(self, pxndot, pxnddt, pxldot, delta):
        _, pxndot, pxnddt, pxldot = self.DeepCalcDotTerms(pxndot, pxnddt, pxldot)

        self._dp_xli    = self._dp_xli + pxldot * delta + pxndot * self._dp_step2
        self._dp_xni    = self._dp_xni + pxndot * delta + pxnddt * self._dp_step2
        self._dp_atime  = self._dp_atime + delta

        return pxndot, pxnddt, pxldot

    def DeepSecular(self, xmdf, omgadf, xnode, emm, xincc, xnn, tsince):
        # Deep space secular effects
        xmdf    = xmdf + self._dp_ssl * tsince
        omgadf  = omgadf + self._dp_ssg * tsince
        xnode   = xnode + self._dp_ssh * tsince
        emm     = self.Orbit.Eccentricity + self._dp_sse * tsince
        xincc   = self.Orbit.Inclination + self._dp_ssi * tsince

        if xincc < 0.0:
            xincc   = -xincc
            xnode   = xnode  + Globals.Pi
            omgadf  = omgadf - Globals.Pi

        xnddt   = 0.0
        xndot   = 0.0
        xldot   = 0.0
        ft      = 0.0
        delt    = 0.0

        fDone   = False

        if self._gp_reso:
            while not fDone:
                if self._dp_atime == 0.0 or \
                    (tsince >= 0.0 and self._dp_atime < 0.0) or \
                    (tsince < 0.0 and self._dp_atime >= 0.0):
                    delt = self._dp_stepn if tsince < 0 else self._dp_stepp

                    # Epoch restart
                    self._dp_atime  = 0.0
                    self._dp_xni    = self._m_xnodp
                    self._dp_xli    = self._dp_xlamo

                    fDone = True
                else:
                    if math.fabs(tsince) < math.fabs(self._dp_atime):
                        delt = self._dp_stepp

                        if tsince >= 0.0:
                            delt = self._dp_stepn

                        xndot, xnddt, xldot = self.DeepCalcIntegrator(xnddt, xnddt, xldot, delt)
                    else:
                        delt = self._dp_stepn

                        if tsince > 0.0:
                            delt = self._dp_stepp

                        fDone = True

            while math.fabs(tsince - self._dp_atime) >= self._dp_stepp:
                xndot, xnddt, xldot = self.DeepCalcIntegrator(xndot, xnddt, xldot, delt)

            ft = tsince - self._dp_atime

            _, xndot, xnddt, xldot = self.DeepCalcDotTerms(xndot, xnddt, xldot)

            xnn  = self._dp_xni + xndot * ft + xnddt * ft * ft * 0.5

            xl   = self._dp_xli + xldot * ft + xndot * ft * ft * 0.5
            temp = -xnode + self._dp_thgr + tsince * NoradSDP4.thdt

            xmdf = xl - omgadf + temp

            if not self._gp_sync:
                xmdf = xl + temp + temp

        return True, xmdf, omgadf, xnode, emm, xincc, xnn, tsince

    def DeepPeriodics(self, e, xincc, omgadf, xnode, xmam, tsince):
        # Lunar-solar periodics
        sinis   = math.sin(xincc)
        cosis   = math.cos(xincc)

        sghs    = 0.0
        shs     = 0.0
        shl     = 0.0
        pe      = 0.0
        pinc    = 0.0
        pl      = 0.0
        sghl    = 0.0

        zm      = self._dp_zmos + NoradSDP4.zns * tsince
        zf      = zm + 2.0 * NoradSDP4.zes * math.sin(zm)
        sinzf   = math.sin(zf)
        f2      =  0.5 * sinzf * sinzf - 0.25
        f3      = -0.5 * sinzf * math.cos(zf)
        ses     = self._dp_se2 * f2 + self._dp_se3 * f3
        sis     = self._dp_si2 * f2 + self._dp_si3 * f3
        sls     = self._dp_sl2 * f2 + self._dp_sl3 * f3 + self._dp_sl4 * sinzf

        sghs    = self._dp_sgh2 * f2 + self._dp_sgh3 * f3 + self._dp_sgh4 * sinzf
        shs     = self._dp_sh2 * f2 + self._dp_sh3 * f3
        zm      = self._dp_zmol + NoradSDP4.znl * tsince
        zf      = zm + 2.0 * NoradSDP4.zel * math.sin(zm)
        sinzf   = math.sin(zf)
        f2      =  0.5 * sinzf * sinzf - 0.25
        f3      = -0.5 * sinzf * math.cos(zf)

        sel     = self._dp_ee2 * f2 + self._dp_e3 * f3
        sil     = self._dp_xi2 * f2 + self._dp_xi3 * f3
        sll     = self._dp_xl2 * f2 + self._dp_xl3 * f3 + self._dp_xl4 * sinzf

        sghl    = self._dp_xgh2 * f2 + self._dp_xgh3 * f3 + self._dp_xgh4 * sinzf
        sh1     = self._dp_xh2 * f2 + self._dp_xh3 * f3
        pe      = ses + sel
        pinc    = sis + sil
        pl      = sls + sll

        pgh     = sghs + sghl
        ph      = shs + shl

        xincc   = xincc + pinc
        e       = e + pe

        if self._dp_xqncl >= 0.2:
            # Apply periodics directly
            ph      = ph / self._m_sinio
            pgh     = pgh - self._m_cosio * ph
            omgadf  = omgadf + pgh
            xnode   = xnode + ph
            xmam    = xmam + pl
        else:
            # Apply periodics with Lyddane modification
            sinok   = math.sin(xnode)
            cosok   = math.cos(xnode)
            alfdp   = sinis * sinok
            betdp   = sinis * cosok
            dalf    =  ph * cosok + pinc * cosis * sinok
            dbet    = -ph * sinok + pinc * cosis * cosok

            alfdp   = alfdp + dalf
            betdp   = betdp + dbet

            xls     = xmam + omgadf + cosis * xnode
            dls     = pl + pgh - pinc * xnode * sinis

            xls     = xls + dls
            xnode   = Globals.AcTan(alfdp, betdp)
            xmam    = xmam + pl
            omgadf  = xls - xmam - math.cos(xincc) * xnode

        return True, e, xincc, omgadf, xnode, xmam

    ##
    # @brief Calculate satellite ECI position/velocity for a given time
    #        This procedure returns the ECI position and velocity for the satellite
    #        in the orbit at the given number of minutes since the TLE epoch time.
    #        The algorithm uses NORAD`s Simplified General Perturbation 4 deep space
    #        orbit model
    #
    # @param tsince Target time, in minutes-past-epoch format. 
    #
    # @return AU-based position/velocity ECI coordinates.
    def GetPosition(self, tsince):
        # Update for secular gravity and atmospheric drag
        xmdf    = self.Orbit.MeanAnomaly + self._m_xmdot * tsince
        omgadf  = self.Orbit.ArgPerigee + self._m_omgdot * tsince
        xnoddf  = self.Orbit.RAAN + self._m_xnodot * tsince
        tsq     = tsince * tsince
        xnode   = xnoddf * self._m_xnodcf * tsq
        tempa   = 1.0 - self._m_c1 * tsince
        tempe   = self.Orbit.BStar * self._m_c4 * tsince
        templ   = self._m_t2cof * tsq
        xn      = self._m_xnodp
        em      = 0.0
        xinc    = 0.0

        _, xmdf, omgadf, xnode, em, xinc, xn, tsince = self.DeepSecular(xmdf, omgadf, xnode, em, xinc, xn, tsince)

        a       = math.pow(Globals.Xke / xn, 2.0 / 3.0) * Globals.Sqr(tempa)
        e       = em - tempe
        xmam    = xmdf + self._m_xnodp * templ

        _, e, xinc, omgadf, xnode, xmam = self.DeepPeriodics(e, xinc, omgadf, xnode, xmam, tsince)

        xl      = xmam + omgadf + xnode

        xn      = Globals.Xke / math.pow(a, 1.5)

        return self.FinalPosition(xinc, omgadf, e, a, xl, xnode, xn, tsince)
