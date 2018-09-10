##
# @file Satellite.py
# @brief 
# @author df_justforfun@163.com
# @version 1.0
# @date 2018-07-18


from pythonOrbitTools.Orbit.Orbit import Orbit

##
# @brief Class to encapsulate a satellite
class Satellite(object):

    # region Properties

    @property
    ##
    # @brief The satellite name.
    #
    # @return 
    def Name(self):
        return self._name

    @property
    ##
    # @brief Information related to the satellite`s orbit.
    #
    # @return 
    def Orbit(self):
        return self._orbit

    # endregion

    ##
    # @brief Standard constructor.
    #
    # @param tle Tle data.
    # @param name Optional satellite name.
    #
    # @return 
    def __init__(self, tle, name=""):
        self._orbit =Orbit(tle)

        if name == "":
            self._name = self.Orbit.SatName
        else:
            self._name = name


    ##
    # @brief Returns the ECI position of the satellite.
    #
    # @param utc The time (UTC) of position calculation.
    #
    # @return The ECI location of the satellite at the given time.
    def PositionEciByDateTime(self, utc):
        return self.Orbit.PositionEciByDateTime(utc)

    ##
    # @brief Returns the ECI position of the satellite.
    #
    # @param mpe The time of position calculation, in minutes-past-epoch.
    #
    # @return The ECI location of the satellite at the given time.
    def PositionEciByMpe(self, mpe):
        return self.Orbit.PositionEciByMpe(mpe)
