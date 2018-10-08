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

    def __init__(self, theta, alpha, a, d, i):
        self.theta = theta
        self.alpha = alpha
        self.d = d
        self.a = a
        self.mass = sp.Symbol('m'+i)
        self.Ixx = sp.Symbol('I' + i + 'xx')
        self.Iyy = sp.Symbol('I' + i + 'yy')
        self.Izz = sp.Symbol('I' + i + 'zz')
        self.intertial_tensor = sp.Matrix([self.Ixx, 0, 0], [0, self.Iyy, 0], [0, 0, self.Izz])
        self.mass_center = sp.Matrix([0, 0, 0])

        if type(type(theta)) == sp.function.UndefinedFunction:
            type = LinkType.rotational
        elif type(type(d)) == sp.function.UndefinedFunction:
            type = LinkType.linear
        else:
            type = LinkType.none
