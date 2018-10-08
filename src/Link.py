import sympy as sp
from enum import Enum


class LinkType(Enum):
    none = -1
    rotational = 0
    linear = 1


class Link:
    type = LinkType.none
    theta = 0
    alpha = 0
    a = 0
    d = 0
    mass = 0
    inertial_tensor = 0

    def __init__(self, theta, alpha, a, d, i, mass_center=None):
        self.theta = theta
        self.alpha = alpha
        self.d = d
        self.a = a
        self.mass = sp.Symbol('m' + str(i))
        self.Ixx = sp.Symbol('I' + str(i) + 'xx')
        self.Iyy = sp.Symbol('I' + str(i) + 'yy')
        self.Izz = sp.Symbol('I' + str(i) + 'zz')
        self.inertial_tensor = sp.Matrix([[self.Ixx, 0, 0], [0, self.Iyy, 0], [0, 0, self.Izz]])
        if type(mass_center) == sp.Matrix:
            self.mass_center = mass_center
        else:
            self.mass_center = sp.Matrix([0, 0, 0])

        if type(type(theta)) == sp.function.UndefinedFunction:
            self.type = LinkType.rotational
        elif type(type(d)) == sp.function.UndefinedFunction:
            self.type = LinkType.linear
        else:
            self.type = LinkType.none
