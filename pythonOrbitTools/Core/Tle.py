##
# @file Tle.py
# @brief 
# @author df_justforfun@163.com
# @version 1.0
# @date 2018-07-16

#########################################################################
#
# NASA Two-Line Element Data format
#
# [Reference: Dr. T.S. Kelso / www.celestrak.com]
#
# Two-line element data consists of three lines in the following format:
#
# AAAAAAAAAAAAAAAAAAAAAAAA
# 1 NNNNNU NNNNNAAA NNNNN.NNNNNNNN +.NNNNNNNN +NNNNN-N +NNNNN-N N NNNNN
# 2 NNNNN NNN.NNNN NNN.NNNN NNNNNNN NNN.NNNN NNN.NNNN NN.NNNNNNNNNNNNNN
#
# Line 0 is a twenty-four-character name.
#
# Lines 1 and 2 are the standard Two-Line Orbital Element Set Format identical
# to that used by NORAD and NASA.  The format description is:
#
#     Line 1
#     Column    Description
#     01-01     Line Number of Element Data
#     03-07     Satellite Number
#     10-11     International Designator (Last two digits of launch year)
#     12-14     International Designator (Launch number of the year)
#     15-17     International Designator (Piece of launch)
#     19-20     Epoch Year (Last two digits of year)
#     21-32     Epoch (Julian Day and fractional portion of the day)
#     34-43     First Time Derivative of the Mean Motion
#               or Ballistic Coefficient (Depending on ephemeris type)
#     45-52     Second Time Derivative of Mean Motion (decimal point assumed;
#               blank if N/A)
#     54-61     BSTAR drag term if GP4 general perturbation theory was used.
#               Otherwise, radiation pressure coefficient.  (Decimal point assumed)
#     63-63     Ephemeris type
#     65-68     Element number
#     69-69     Check Sum (Modulo 10)
#               (Letters, blanks, periods, plus signs = 0; minus signs = 1)
#
#     Line 2
#     Column    Description
#     01-01     Line Number of Element Data
#     03-07     Satellite Number
#     09-16     Inclination [Degrees]
#     18-25     Right Ascension of the Ascending Node [Degrees]
#     27-33     Eccentricity (decimal point assumed)
#     35-42     Argument of Perigee [Degrees]
#     44-51     Mean Anomaly [Degrees]
#     53-63     Mean Motion [Revs per day]
#     64-68     Revolution number at epoch [Revs]
#     69-69     Check Sum (Modulo 10)
#
#     All other columns are blank or fixed.
#
# Example:
#
# NOAA 6
# 1 11416U          86 50.28438588 0.00000140           67960-4 0  5293
# 2 11416  98.5105  69.3305 0012788  63.2828 296.9658 14.24899292346978

from enum import Enum, unique
from pythonOrbitTools.Core.Globals import Globals
from pythonOrbitTools.Core.Julian import Julian

@unique
class Line(Enum):
    Zero = 0    # Satellite Name
    One  = 1    # The first line of TLE
    Two  = 2    # The Second line of TLE

@unique
class Field(Enum):
    NoradNum       = 0
    IntlDesc       = 1
    SetNumber      = 2     # TLE set number
    EpochYear      = 3     # Epoch: Last two digits of year
    EpochDay       = 4     # Epoch: Fractional Julian Day of year
    OrbitAtEpoch   = 5     # Orbit at epoch
    Inclination    = 6     # Inclination
    Raan           = 7     # R.A. ascending node
    Eccentricity   = 8     # Eccentricity
    ArgPerigee     = 9     # Argument of perigee
    MeanAnomaly    = 10    # Mean anomaly
    MeanMotion     = 11    # Mean motion
    MeanMotionDt   = 12    # First time derivative of mean motion
    MeanMotionDt2  = 13    # Second time derivative of mean motion
    BStarDrag      = 14    # BSTAR Drag

@unique
class Unit(Enum):
    Radians = 0     # radians
    Degrees = 1     # degrees
    Native  = 2     # TLE format native units (no conversion)    



##
# @brief This class encapsulates a single set of standard NORAD two-line elements.
class Tle(object):

    #region Column Offsets

    # Note: The column offsets are zero-based.

    # Name
    TLE_LEN_LINE_DATA      = 69; TLE_LEN_LINE_NAME      = 24;

    # Line 1
    TLE1_COL_SATNUM        =  2; TLE1_LEN_SATNUM        =  5
    TLE1_COL_INTLDESC_A    =  9; TLE1_LEN_INTLDESC_A    =  2
    TLE1_COL_INTLDESC_B    = 11; TLE1_LEN_INTLDESC_B    =  3
    TLE1_COL_INTLDESC_C    = 14; TLE1_LEN_INTLDESC_C    =  3
    TLE1_COL_EPOCH_A       = 18; TLE1_LEN_EPOCH_A       =  2
    TLE1_COL_EPOCH_B       = 20; TLE1_LEN_EPOCH_B       = 12
    TLE1_COL_MEANMOTIONDT  = 33; TLE1_LEN_MEANMOTIONDT  = 10
    TLE1_COL_MEANMOTIONDT2 = 44; TLE1_LEN_MEANMOTIONDT2 =  8
    TLE1_COL_BSTAR         = 53; TLE1_LEN_BSTAR         =  8
    TLE1_COL_EPHEMTYPE     = 62; TLE1_LEN_EPHEMTYPE     =  1
    TLE1_COL_ELNUM         = 64; TLE1_LEN_ELNUM         =  4

    # Line 2
    TLE2_COL_SATNUM        = 2;  TLE2_LEN_SATNUM        =  5
    TLE2_COL_INCLINATION   = 8;  TLE2_LEN_INCLINATION   =  8
    TLE2_COL_RAASCENDNODE  = 17; TLE2_LEN_RAASCENDNODE  =  8
    TLE2_COL_ECCENTRICITY  = 26; TLE2_LEN_ECCENTRICITY  =  7
    TLE2_COL_ARGPERIGEE    = 34; TLE2_LEN_ARGPERIGEE    =  8
    TLE2_COL_MEANANOMALY   = 43; TLE2_LEN_MEANANOMALY   =  8
    TLE2_COL_MEANMOTION    = 52; TLE2_LEN_MEANMOTION    = 11
    TLE2_COL_REVATEPOCH    = 63; TLE2_LEN_REVATEPOCH    =  5

    #endregion


    ##
    # @brief 
    #
    # @param strLine1: the first line of TLE
    # @param strLine2: the second line of TLE
    # @param strName: satellite name, default value is ""
    #
    # @return 
    def __init__(self, strLine1, strLine2, strName=""):
        self._m_Line0 = strName
        self._m_Line1 = strLine1
        self._m_Line2 = strLine2

        # Converted fields, in float()-able form
        # Key   - Filed
        # Value - String
        self._m_Field = {}
        # Cache of field values in "float" format
        # Key   - integer
        # Value - float
        self._m_Cache = {}

        self.Initialize()

    # region Properties

    @property
    def Name(self):
        return self._m_Line0

    @property
    def Line1(self):
        return self._m_Line1

    @property
    def Line2(self):
        return self._m_Line2

    @property
    def NoradNum(self):
        return self.GetFieldAsString(Field.NoradNum, False)

    @property
    def Eccentricity(self):
        return self.GetFieldAsString(Field.Eccentricity, False)

    @property
    def Inclination(self):
        return self.GetFieldAsString(Field.Inclination, True)

    @property
    def Epoch(self):
        return "{0:02.0f}{1:012.8f}".format(self.GetFieldAsValue(Field.EpochYear),
                                            self.GetFieldAsValue(Field.EpochDay))

    @property
    def IntlDesciption(self):
        return self.GetFieldAsString(Field.IntlDesc, False)

    @property
    def SetNumber(self):
        return self.GetFieldAsString(Field.SetNumber, False)

    @property
    def OrbitAtEpoch(self):
        return self.GetFieldAsString(Field.OrbitAtEpoch, False)

    @property
    def RAAscendingNode(self):
        return self.GetFieldAsString(Field.Raan, True)

    @property
    def ArgPerigee(self):
        return self.GetFieldAsString(Field.ArgPerigee, True)

    @property
    def MeanAnomaly(self):
        return self.GetFieldAsString(Field.MeanAnomaly, True)

    @property
    def MeanMotion(self):
        return self.GetFieldAsString(Field.MeanMotion, True)

    @property
    def MeanMotionDt(self):
        return self.GetFieldAsString(Field.MeanMotionDt, False)

    @property
    def MeanMotionDt2(self):
        return self.GetFieldAsString(Field.MeanMotionDt2, False)

    @property
    def BStarDrag(self):
        return self.GetFieldAsString(Field.BStarDrag, False)

    @property
    def EpochJulian(self):
        epochYear = int(self.GetFieldAsValue(Field.EpochYear))
        epochDay  = self.GetFieldAsValue(Field.EpochDay)

        if epochYear < 57:
            epochYear += 2000
        else:
            epochYear += 1900

        return Julian().InitializeByYearAndDoy(epochYear, epochDay)

    # end region

    def Initialize(self):
        self._m_Field[Field.NoradNum] = self._m_Line1[Tle.TLE1_COL_SATNUM:
                                                      Tle.TLE1_COL_SATNUM+Tle.TLE1_LEN_SATNUM]
        self._m_Field[Field.IntlDesc] = self._m_Line1[Tle.TLE1_COL_INTLDESC_A:
                                                      Tle.TLE1_COL_INTLDESC_A+
                                                      Tle.TLE1_LEN_INTLDESC_A+
                                                      Tle.TLE1_LEN_INTLDESC_B+
                                                      Tle.TLE1_LEN_INTLDESC_C]
        self._m_Field[Field.EpochYear] = self._m_Line1[Tle.TLE1_COL_EPOCH_A:
                                                       Tle.TLE1_COL_EPOCH_A+Tle.TLE1_LEN_EPOCH_A]
        self._m_Field[Field.EpochDay] = self._m_Line1[Tle.TLE1_COL_EPOCH_B:
                                                      Tle.TLE1_COL_EPOCH_B+Tle.TLE1_LEN_EPOCH_B]

        if self._m_Line1[Tle.TLE1_COL_MEANMOTIONDT] == "-":
            # value is negative
            self._m_Field[Field.MeanMotionDt] = "-0"
        else:
            self._m_Field[Field.MeanMotionDt] = "0"

        self._m_Field[Field.MeanMotionDt] += self._m_Line1[Tle.TLE1_COL_MEANMOTIONDT+1:
                                                           Tle.TLE1_COL_MEANMOTIONDT+1+Tle.TLE1_LEN_MEANMOTIONDT]

        # decimal point assumed; exponential notation
        self._m_Field[Field.MeanMotionDt2] = Tle.ExpToDecimal(self._m_Line1[Tle.TLE1_COL_MEANMOTIONDT2:
                                                                            Tle.TLE1_COL_MEANMOTIONDT2+Tle.TLE1_LEN_MEANMOTIONDT2])

        # decimal point assumed; exponential notaion
        self._m_Field[Field.BStarDrag] = Tle.ExpToDecimal(self._m_Line1[Tle.TLE1_COL_BSTAR:
                                                                        Tle.TLE1_COL_BSTAR+Tle.TLE1_LEN_BSTAR])

        # TLE1_COL_EPHEMTYPE
        # TLE1_LEN_EPHEMTYPE

        self._m_Field[Field.SetNumber] = self._m_Line1[Tle.TLE1_COL_ELNUM:
                                                       Tle.TLE1_COL_ELNUM+Tle.TLE1_LEN_ELNUM].lstrip()

        # TLE2_COL_SATNUM
        # TLE2_LEN_SATNUM

        self._m_Field[Field.Inclination] = self._m_Line2[Tle.TLE2_COL_INCLINATION:
                                                         Tle.TLE2_COL_INCLINATION+Tle.TLE2_LEN_INCLINATION].lstrip()
        self._m_Field[Field.Raan] = self._m_Line2[Tle.TLE2_COL_RAASCENDNODE:
                                                  Tle.TLE2_COL_RAASCENDNODE+Tle.TLE2_LEN_RAASCENDNODE].lstrip()
        # Eccentricity: decimal point is assumed
        self._m_Field[Field.Eccentricity] = "0."+self._m_Line2[Tle.TLE2_COL_ECCENTRICITY:
                                                               Tle.TLE2_COL_ECCENTRICITY+Tle.TLE2_LEN_ECCENTRICITY]
        self._m_Field[Field.ArgPerigee] = self._m_Line2[Tle.TLE2_COL_ARGPERIGEE:
                                                        Tle.TLE2_COL_ARGPERIGEE+Tle.TLE2_LEN_ARGPERIGEE].lstrip()
        self._m_Field[Field.MeanAnomaly] = self._m_Line2[Tle.TLE2_COL_MEANANOMALY:
                                                         Tle.TLE2_COL_MEANANOMALY+Tle.TLE2_LEN_MEANANOMALY].lstrip()
        self._m_Field[Field.MeanMotion] = self._m_Line2[Tle.TLE2_COL_MEANMOTION:
                                                        Tle.TLE2_COL_MEANMOTION+Tle.TLE2_LEN_MEANMOTION].lstrip()
        self._m_Field[Field.OrbitAtEpoch] = self._m_Line2[Tle.TLE2_COL_REVATEPOCH:
                                                          Tle.TLE2_COL_REVATEPOCH+Tle.TLE2_LEN_REVATEPOCH].lstrip()


    ##
    # @brief: Returns the requested TLE data field as a type double. 
    #         The numeric return values are cached; requesting the same field repeatedly incurs minimal overhead.
    #
    # @param fld: The TLE field to retrievv.
    # @param units: Specifies the units desired.
    #
    # @return: The requested field`s value, converted to the correct units if necessary. 
    def GetFieldAsValue(self, fld, units=Unit.Native):
        # Return cache contents if it exists, else populate cache.
        key = Tle.Key(units, fld)

        if key in self._m_Cache:
            # return cached value
            return float(self._m_Cache[key])
        else:
            # Value not in cache; add it
            valNative           = float(self._m_Field[fld])
            valConv             = Tle.ConvertUnits(valNative, fld, units)
            self._m_Cache[key]  = valConv

            return valConv


    ##
    # @brief Returns the requested TLE data field in native form as a text string
    #
    # @param fld: The TLE field to retrieve.
    # @param appendUnits: If true, the native units are appended to the end of the returned string.
    #
    # @return The requested field as a string.
    def GetFieldAsString(self, fld, appendUnits):

        string = self._m_Field[fld]

        if appendUnits:
            string += Tle.GetUnits(fld)

        return string.strip()


    @staticmethod
    ##
    # @brief Generates a key for the TLE field cache.
    #
    # @param u: The data type is Unit(Enum)
    # @param f: The data type is Field(Enum)
    #
    # @return 
    def Key(u, f):
        return u.value * 100 + f.value


    @staticmethod
    ##
    # @brief: Converts the given TLE field to the request units. 
    #
    # @param valNative: Value to convert (native units).
    # @param fld: Field ID of the value being converted.
    # @param units: Units to convert to.
    #
    # @return: The converted value. 
    def ConvertUnits(valNative, fld, units):
        if fld == Field.Inclination or fld == Field.Raan or fld == Field.ArgPerigee or fld == Field.MeanAnomaly:
            # The native TLE format is degrees
            if units == Unit.Radians:
                return Globals.ToRadians(valNative)

        return valNative



    @staticmethod
    def GetUnits(fld):
        strDegrees      = " degrees"
        strRevsPerDay   = " revs / day"

        if fld == Field.Inclination or fld == Field.Raan or fld == Field.ArgPerigee or fld == Field.MeanAnomaly:
               return strDegrees

        if fld == Field.MeanMotion:
            return strRevsPerDay
        else:
            return ""



    ##
    # @brief 
    # Converts TLE-style exponential notation of the form [ |+|-]00000[ |+|-]0
    # to decimal notaion, Assumes implied decimal point to the left of the first
    # number in the string, i.e.,
    #       " 12345-3" =  0.00012345
    #       "-23429-5" = -0.0000023429 
    #       " 40436+1" = 4.0436
    # Also assumes that lack of a sign character implies a positive value, i.e.,
    #       " 00000 0" =  0.00000
    #       " 31415 1" =  3.1415
    @staticmethod
    def ExpToDecimal(string):

        COL_SIGN     = 0
        LEN_SIGN     = 1

        COL_MANTISSA = 1
        LEN_MANTISSA = 5

        COL_EXPONENT = 6
        LEN_EXPONENT = 2

        sign        = string[COL_SIGN:      COL_SIGN+LEN_SIGN]
        mantissa    = string[COL_MANTISSA:  COL_MANTISSA+LEN_MANTISSA]
        exponent    = string[COL_EXPONENT:  COL_EXPONENT+LEN_EXPONENT].lstrip()

        val = float(sign+"0."+mantissa+"e"+exponent)

        return str(val)

    @staticmethod
    ##
    # @brief: Determines if a given string has the expected format of a single 
    #         line of TLE data.
    #         A valid satellite name is less than or equal to TLE_LEN_LINE_NAME
    #           characters;
    #         A valid data line must:
    #           Have as the first character the line number
    #           Have as the second character a blank
    #           Be TLE_LEN_LINE_DATA characters long
    #           etc.
    #
    # @param string: The input string.
    # @param line: The line ID of the input string. (Line.Zero, Line.One, Line.Two)
    #
    # @return: True if the input string has the format of the given line ID. 
    def IsValidFormat(string, line):
        raise NotImplementedError


    @staticmethod
    ##
    # @brief: Calculate the check sum for a given line of TLE data, the last character 
    #         of which is the current checksum. (Although there is no check here,
    #         the current checksum should match the one calculated.)
    #         The checksum algorithm:
    #           Each number in the data line is summed, module 10.
    #           Non-numeric characters are zero, except minus signs, which are 1.
    #
    # @param string: The input string.
    #
    # @return 
    def CheckSum(string):
        raise NotImplementedError

if __name__ == "__main__":

    line1 = "1 25994U 99068A   18196.75093423 -.00000025  00000-0  45345-5 0  9993"
    line2 = "2 25994  98.2051 271.2050 0001021  68.8940 291.2371 14.57112414987988"

    tle = Tle(line1, line2, "TERRA")

    print("SateName             : {}".format(tle.Name))
    print("Line1                : {}".format(tle.Line1))
    print("Line2                : {}".format(tle.Line2))
    print("NoradNumber          : {}".format(tle.NoradNum))
    print("Eccentricity(str)    : {}".format(tle.Eccentricity))
    print("Eccentricity(val)    : {}".format(tle.GetFieldAsValue(Field.Eccentricity)))
    print("Inclination(str)     : {}".format(tle.Inclination))
    print("Inclination(val)     : {} rad".format(tle.GetFieldAsValue(Field.Inclination, Unit.Radians)))
    print("Epoch                : {}".format(tle.Epoch))
    print("IntlDesciption       : {}".format(tle.IntlDesciption))
    print("SetNumber            : {}".format(tle.SetNumber))
    print("OrbitAtEpoch         : {}".format(tle.OrbitAtEpoch))
    print("RAAscendingNode(str) : {}".format(tle.RAAscendingNode))
    print("RAAscendingNode(val) : {} rad".format(tle.GetFieldAsValue(Field.Raan, Unit.Radians)))
    print("ArgPerigee(str)      : {}".format(tle.ArgPerigee))
    print("ArgPerigee(val)      : {} rad".format(tle.GetFieldAsValue(Field.ArgPerigee, Unit.Radians)))
    print("MeanAnomaly(str)     : {}".format(tle.MeanAnomaly))
    print("MeanAnomaly(val)     : {} rad".format(tle.GetFieldAsValue(Field.MeanAnomaly, Unit.Radians)))
    print("MeanMotion(str)      : {}".format(tle.MeanMotion))
    print("MeanMotion(val)      : {}".format(tle.GetFieldAsValue(Field.MeanMotion)))
    print("MeanMotionDt(str)    : {}".format(tle.MeanMotionDt))
    print("MeanMotionDt(val)    : {}".format(tle.GetFieldAsValue(Field.MeanMotionDt)))
    print("MeanMotionDt2        : {}".format(tle.MeanMotionDt2))
    print("BStarDrag(str)       : {}".format(tle.BStarDrag))
    print("BStarDrag(val)       : {}".format(tle.GetFieldAsValue(Field.BStarDrag)))
    print("EpochJulian          : {}".format(tle.EpochJulian.Date))
    print("EpochYear            : {}".format(tle.GetFieldAsValue(Field.EpochYear)))
    print("EpochDay             : {}".format(tle.GetFieldAsValue(Field.EpochDay)))
    print("EpochDateTime        : {}".format(tle.EpochJulian.ToTime()))
