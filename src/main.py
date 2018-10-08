import math
from robotlib import *
import sympy as sp
from Robot import Robot

# definicja symboli
# definicja linków
# dodanie linków do robota
# kalkulacja dynamiki

t = sp.Symbol('t')

q1 = sp.Function('q1')(t)
q2 = sp.Function('q2')(t)
q3 = sp.Function('q3')(t)

a1 = sp.Function('a1')(t)
a2 = sp.Function('a2')(t)

c1 = sp.Matrix([0, 0, sp.Symbol('c1z')])
c2 = sp.Matrix([0, 0, sp.Symbol('c2z')])
c3 = sp.Matrix([0, 0, sp.Symbol('c3z')])
# theta, alpha , a , d
link1 = Link(q1, 0.0, 0.0, 0.0, 1, c1)
link2 = Link(q2, -math.pi/2.0, 0.0, 0.0, 2, c2)
link3 = Link(0.0, 0.0, q1, 0.0, 3, c3)

my_robot = Robot(t)
my_robot.add_link(link1)
my_robot.add_link(link2)
my_robot.add_link(link3)

my_robot.calculate_dynamics()

my_robot.print_transofrmation_matrices()
my_robot.print_transofrmation_matrices_in_base()

my_robot.print_veocitioes()

my_robot.print_dynamics()

M = my_robot.M

sp.pprint(M[0,0].diff(q2))
