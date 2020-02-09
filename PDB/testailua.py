from sys import argv
import os
from math import sqrt, atan2, degrees

def dot_product(v1, v2):
    """ Calculate the dot product of two vectors """
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]

def cross_product(v1, v2):
    """ Calculate the cross product of two vectors """
    i = v1[1] * v2[2] - v1[2] * v2[1]
    j = v1[2] * v2[0] - v1[0] * v2[2]
    k = v1[0] * v2[1] - v1[1] * v2[0]
    return [i, j, k]

def magnitude(v):
    """ Calculate the size of a vector """
    return sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2)

def calculateDihedral(a1, a2, a3, a4):
    """ Calculates the normal vector of the planes
    defined by four atom coordinates """
    # START CODING HERE
    # calculate normal vectors to the planes defined by a1,a2,a3 and a2,a3,a4
    # you may use the functions "cross_product","dot_product" and "magnitude" defined above
    # you can also use the python math function "atan2" and "degrees"

    # Calculate vectors and normals from the points
    #a1 = cn
    #a2 = nca
    #a3 = nca
    #a4 = cac
    
    n1 = cross_product(a2,a1)
    n2 = cross_product(a3,a1)
    n3 = cross_product(a3,a2)
    n4 = cross_product(a4,a2)

    #cos
    sin = dot_product(n1,n2) / (magnitude(n1) * magnitude(n2))
    cos = dot_product(n3,n4) / (magnitude(n3) * magnitude(n4))

    #atan2
    dihedral = degrees(atan2(sin, cos))

    # END CODING HERE
    return dihedral

print(calculateDihedral([1, 9, 2], [3, 2, 1], [2, 4, 7], [8, 2, 5]))