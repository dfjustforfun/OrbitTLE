##
# @file Julian.py
# @brief 
#
#  This class encapsulates a Julian date system where the day starts at noon.
#  Some Julian dates:
#       01/01/1990 00:00 UTC - 2447892.5
#       01/01/1990 12:00 UTC - 2447893.0
#       01/01/2000 00:00 UTC - 2451544.5
#       01/01/2001 00:00 UTC - 2451910.5
#  The Julian day begins at noon, which allows astronomeers to have the
#  same date in a single observing session.
#  
#  References:
#       "Astronomical Formulate for Calculators", Jean Meenus, 4th Edition     
#       "Satellite Communications", Dennis Roddy, 2nd Edition, 1995
#       "Spacecraft Attitude Determination and Control", James R. Wertz, 1984
#
# @author df_justforfun@163.com
# @version 1.0
# @date 2018-07-18

import datetime
from pythonOrbitTools.Core.Globals import Globals

##
# @brief Encapsulates a Julian date.
class Julian(object):

    EPOCH_JAN0_12H_1900 = 2415020.0  # Dec 31.5 1899 = Dec 31 1899 12h UTC
    EPOCH_JAN1_00H_1900 = 2415020.5  # Jan  1.0 1900 = Jan  1 1900 00h UTC
    EPOCH_JAN1_12H_1900 = 2415021.0  # Jan  1.5 1900 = Jan  1 1900 12h UTC
    EPOCH_JAN1_12H_2000 = 2451545.0  # Jan  1.5 2000 = Jan  1 2000 12h UTC


    def __init__(self):

        self._m_Date    = 0.0 # Julian_date
        self._m_Year    = 0   # Year including century
        self._m_Day     = 0.0 # Day of year, 1.0 = Jan 1 00h


    # region Properties

    @property
    def Date(self):
        return self._m_Date

    @property
    def FromJan0_12h_1900(self):
        return self._m_Date - Julian.EPOCH_JAN0_12H_1900

    @property
    def FromJan1_00h_1900(self):
        return self._m_Date - Julian.EPOCH_JAN1_00H_1900

    @property
    def FromJan1_12h_1900(self):
        return self._m_Date - Julian.EPOCH_JAN1_12H_1900

    @property
    def FromJan1_12h_2000(self):
        return self._m_Date - Julian.EPOCH_JAN1_12H_2000

    # endregion


    ##
    # @brief Initialize the Julian date object.
    #
    # @param utc The UTC time to convert.
    #
    # @return 
    def InitializeByUTC(self, utc):

        day = utc.timetuple().tm_yday + (utc.hour + (utc.minute + (utc.second + (utc.microsecond/1000000.0))/60.0)/60.0)/24.0
        self.InitializeByYearAndDoy(utc.year, day)

        return self


    ##
    # @brief Initialize the Julian date object.
    #        The first day of the year, Jan 1, is day 1.0. Noon on Jan 1 is
    #        represented by the day value of 1.5, etc
    #        The fraction part of the day value is the fraction portion of
    #        the day.
    #        Examples:
    #           day = 1.0  Jan 1 00h
    #           day = 1.5  Jan 1 12h
    #           day = 2.0  Jan 2 00h
    #
    # @param year The year, including the century (i.e., 2012).
    # @param doy Day of year (1 means January 1, etc.)
    #
    # @return 
    def InitializeByYearAndDoy(self, year, doy):

        # Arbitrary years used for error checking
        if year < 1900 or year > 2100:
            raise ValueError("year")

        # The last day of a leap year is day 366
        if doy < 1.0 or doy >= 367.0:
            raise ValueError("doy")

        self._m_Year = year
        self._m_Day  = doy

        # Now calculate Julian date
        # Ref: "Astronomical Formulate for Calculators", Jean Meeus, pages 23-25

        year -= 1

        # Centuries are not leap years unless they divide by 400
        A = int(year/100)
        B = 2 - A + int(A/4)

        NewYears = int(365.25 * year) + int(30.6001 * 14) + 1720994.5 + B

        self._m_Date = NewYears + doy

        return self

    ##
    # @brief Calculates the time difference between two Julian dates.
    #
    # @param date Julian date.
    #
    # @return A timespan representing the time difference between the two dates.
    def Diff(self, date):
        return datetime.timedelta(days=(self._m_Date - date.Date))


    ##
    # @brief Calculate Greenwich Mean Sidereal Time for the Julian date.
    #
    # @return The angle, in radians, measuring eastward from the Vernal Equinox to 
    #         the prime meridian. This angle is also referred to as "ThetaG"(Theta GMST)
    def ToGmst(self):
        # References:
        #   The 1992 Astronomical Almanac, page B6.
        #   Explanatory Supplement to the Astronomical Almanac, page 50.
        #   Orbital Coordinate Systems, Part III, Dr. T.S. Kelso,
        #       Satellite Times, Nov/Dec 1995

        UT = (self._m_Date + 0.5) % 1
        TU = (self.FromJan1_12h_2000-UT)/36525.0

        GMST = 24110.54841 + TU * (8640184.812866 + TU * (0.093104- TU*6.2e-06))

        GMST = (GMST + Globals.SecPerDay * Globals.OmegaE * UT) % Globals.SecPerDay

        if GMST < 0.0:
            GMST += Globals.SecPerDay # "wrap" negative modulo value

        return (Globals.TwoPi * (GMST / Globals.SecPerDay))

    ##
    # @brief Calculate Local Mean Sidereal Time for this Julian date at the given longitude. 
    #
    # @param lon The longitude, in radians, measured west from Greenwich.
    #
    # @return The angle, in radians, measuring eastward from the Vernal Equinox to the given longitude. 
    def ToLmst(self, lon):
        return (self.ToGmst()+lon) % Globals.TwoPi

    ##
    # @brief Returns a UTC DateTime object that corresponds to this Julian date.
    #
    # @return A DateTime object in UTC. 
    def ToTime(self):

        # Jan 1
        dt = datetime.datetime(self._m_Year, 1, 1)

        # _m_Day = 1 = Jan1
        dt = dt + datetime.timedelta(days=(self._m_Day - 1.0))

        return dt

if __name__ == "__main__":
    utcnow = datetime.datetime.utcnow()
    print("utcnow           : {}".format(utcnow))
    julian = Julian().InitializeByUTC(utcnow)
    print("Julian.ToTime()  : {}".format(julian.ToTime()))
    print("Gmst             : {}".format(julian.ToGmst()))
