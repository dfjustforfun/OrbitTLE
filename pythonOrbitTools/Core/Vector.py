##
# @file Vector.py
# @brief 
# @author df_justforfun@163.com
# @version 1.0
# @date 2018-07-18

import math

##
# @brief Encapsulates a simple 4-component vector
class Vector(object):



    ##
    # @brief Creates a new vector with the give XYZ-W components.
    #
    # @param x: The X component.
    # @param y: The Y component.
    # @param z: The Z component.
    # @param w: The W component.
    #
    # @return 
    def __init__(self, x, y, z, w=0.0):

        self._x = x
        self._y = y
        self._z = z
        self._w = w

    # region Properties

    @property
    def X(self):
        return self._x

    @property
    def Y(self):
        return self._y

    @property
    def Z(self):
        return self._z

    @property
    def W(self):
        return self._w

    # endregion


    ##
    # @brief: Multiply each component in the vector by a given factor. 
    #
    # @param factor: The factor.
    #
    # @return 
    def Mul(self, factor):
        self._x *= factor
        self._y *= factor
        self._z *= factor
        self._w *= math.fabs(factor)

    ##
    # @brief: Subtracts a vector from this vector .
    #
    # @param vec: The vector to subtract.
    #
    # @return 
    def Sub(self, vec):
        self._x -= vec.X
        self._y -= vec.Y
        self._z -= vec.Z
        self._w -= vec.W


    ##
    # @brief: Calculates the angle, in randians, between this vector and another 
    #
    # @param vec: The second vector.
    #
    # @return: The angle between the two vectors, in radians 
    def Angle(self, vec):
        return math.acos(self.Dot(vec) / (self.Magnitude() * vec.Magnitude()))

    ##
    # @brief: Calculates the magnitude of the vector. 
    #
    # @return: The vector magnitude. 
    def Magnitude(self):
        return math.sqrt((self._x*self._x) + (self._y*self._y)+ (self._z*self._z))

    ##
    # @brief: Calculates the dot product of this vector and another. 
    #
    # @param vec: The second vector
    #
    # @return: The dot product. 
    def Dot(self, vec):
        return (self._x*vec.X) + (self._y*vec.Y) + (self._z*vec.Z)

    ##
    # @brief: Calculates the distance between two vectors as point in XYZ space. 
    #
    # @param vec: The second vector.
    #
    # @return: The calculated distance. 
    def Distance(self, vec):
        return math.sqrt(math.pow(self._x-vec.X, 2.0) + math.pow(self._y-vec.Y, 2.0) + math.pow(self._z-vec.Z, 2.0))


    ##
    # @brief: Rotates the XYZ coordinates around the X-axis 
    #
    # @param radians
    #
    # @return 
    def RotateX(self, radians):
        y = self._y

        self._y = (math.cos(radians)*y) - (math.sin(radians)*self._z)
        self._z = (math.sin(radians)*y) + (math.cos(radians)*self._z)

    ##
    # @brief: Rotates the XYZ coordinates around the Y-axis 
    #
    # @param radians
    #
    # @return 
    def RotateY(self, radians):
        x = self._x

        self._x = ( math.cos(radians)*x) + (math.sin(radians)*self._z)
        self._z = (-math.sin(radians)*x) + (math.cos(radians)*self._z)

    ##
    # @brief: Rotates the XYZ coordinates around the Z-axis. 
    #
    # @param radians
    #
    # @return 
    def RotateZ(self, radians):
        x = self._x

        self._x = (math.cos(radians)*x) - (math.sin(radians)*self._y)
        self._y = (math.sin(radians)*x) + (math.cos(radians)*self._y)

    ##
    # @brief: Offset the XYZ coordinates. 
    #
    # @param 
    #
    # @return 
    def Translate(self, x, y, z):
        self._x += x
        self._y += y
        self._z += z

if __name__ == "__main__":

    print("******************** Mul ***********************")
    v = Vector(1, 2, 3, 4)
    factor = -1.5
    print("({}, {}, {}, {}) * {}".format(v.X, v.Y, v.Z, v.W, factor))
    v.Mul(factor)
    print("result: ({}, {}, {}, {})\n".format(v.X, v.Y, v.Z, v.W))


    print("******************** Sub ***********************")
    v1 = Vector(4, 5, 6, 7)
    v2 = Vector(4, 3, 2, 1)
    print("({}, {}, {}, {}) - ({}, {}, {}, {})".format(v1.X, v1.Y, v1.Z, v1.W, v2.X, v2.Y, v2.Z, v2.W))
    v1.Sub(v2)
    print("result: ({}, {}, {}, {})\n".format(v1.X, v1.Y, v1.Z, v1.W))


    print("******************** Angle ***********************")
    v1 = Vector(3, 2, 1)
    v2 = Vector(2, 1, 3)
    print("({}, {}, {}) and ({}, {}, {})".format(v1.X, v1.Y, v1.Z, v2.X, v2.Y, v2.Z))
    print("result: {} rad\n".format(v1.Angle(v2)))



    print("******************** Magnitude ***********************")
    v = Vector(3, 4, 5)
    print("Magnitude ({}, {}, {})".format(v.X, v.Y, v.Z))
    print("result: {}\n".format(v.Magnitude()))


    print("******************** Dot ***********************")
    v1 = Vector(2, 3, 1)
    v2 = Vector(-2, 1, 4)
    print("({}, {}, {}) dot ({}, {}, {})".format(v1.X, v1.Y, v1.Z, v2.X, v2.Y, v2.Z))
    print("result: {}\n".format(v1.Dot(v2)))


    print("******************** Distance ***********************")
    v1 = Vector(2, 3, 1)
    v2 = Vector(-2, 1, 4)
    print("({}, {}, {}) distance ({}, {}, {})".format(v1.X, v1.Y, v1.Z, v2.X, v2.Y, v2.Z))
    print("result: {}\n".format(v1.Distance(v2)))


    print("******************** RotateX ***********************")
    v = Vector(2, 2, 2)
    radians = math.pi/4
    print("({}, {}, {}) RotateX {} rad ({} degrees)".format(v.X, v.Y, v.Z, radians, math.degrees(radians)))
    v.RotateX(radians)
    print("result: ({}, {}, {})\n".format(v.X, v.Y, v.Z))


    print("******************** RotateY ***********************")
    v = Vector(2, 2, 2)
    radians = math.pi/4
    print("({}, {}, {}) RotateY {} rad ({} degrees)".format(v.X, v.Y, v.Z, radians, math.degrees(radians)))
    v.RotateY(radians)
    print("result: ({}, {}, {})\n".format(v.X, v.Y, v.Z))


    print("******************** RotateZ ***********************")
    v = Vector(2, 2, 2)
    radians = math.pi/4
    print("({}, {}, {}) RotateZ {} rad ({} degrees)".format(v.X, v.Y, v.Z, radians, math.degrees(radians)))
    v.RotateZ(radians)
    print("result: ({}, {}, {})\n".format(v.X, v.Y, v.Z))


    print("******************** Translate ***********************")
    v = Vector(2, 2, 2)
    print("({}, {}, {}) translate (-1, 3, 1.5)".format(v.X, v.Y, v.Z))
    v.Translate(-1, 3, 1.5)
    print("result: ({}, {}, {})\n".format(v.X, v.Y, v.Z))
