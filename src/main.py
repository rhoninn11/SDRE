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

# a2 = sp.Symbol('a2')
a2 = 0.125
# a3 = sp.Symbol('a3')
a3 = 0.248

# d1 = sp.Symbol('d1')
d1 = 0.689

c1 = sp.Matrix([0, 0, -0.35])
c2 = sp.Matrix([0.06, 0, 0])
c3 = sp.Matrix([0.06, 0, 0])
# theta, alpha , a , d
link1 = Link(q1, 0.0, 0.0, d1, 1, c1, 1.0)
link2 = Link(q2, -math.pi / 2.0, 0.0, 0.0, 2, c2, 1.0)
link3 = Link(q3, 0.0, a2, 0.0, 3, c3, 1.0)
link4 = Link(0.0, 0.0, a3, 0.0, 4)

my_robot = Robot(t)
my_robot.add_link(link1)
my_robot.add_link(link2)
my_robot.add_link(link3)
my_robot.add_link(link4)

my_robot.calculate_dynamics()

my_robot.print_transofrmation_matrices()
my_robot.print_transofrmation_matrices_in_base()

my_robot.print_veocitioes()

my_robot.print_dynamics()

# my_robot.print_to_latex()
#
# my_robot.print_M_elements_to_latex()
# my_robot.print_C_elements_to_latex()
# my_robot.print_G_elements_to_latex()

sim = my_robot.setup_simulation()
sim.simulate(sp.Matrix([1.0, 1.0, 1.0]), sp.Matrix([2.0, 2.0, 2.0]))

