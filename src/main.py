import math
from robotlib import *
import sympy as sp
import sympy.printing.latex as spltx

t = sp.Symbol('t');
q1 = sp.Function('q1')(t)
q2 = sp.Function('q2')(t)
q3 = sp.Function('q3')(t)
cords = sp.Matrix([q1, q2, q3])

t_1in0 = calculate_zdh_matrix(q1, 0.0, 0.0, 0.689)
t_2in1 = calculate_zdh_matrix(q2, -math.pi / 2, 0.0, 0.0)
t_3in2 = calculate_zdh_matrix(q3, 0.0, 0.125, 0.0)
t_4in3 = calculate_zdh_matrix(0.0, 0.0, 0.248, 0.0)

t_4in0 = t_1in0 * t_2in1 * t_3in2 * t_4in3
t_3in0 = t_1in0 * t_2in1 * t_3in2
t_2in0 = t_1in0 * t_2in1

# print('T 1 in 0-------------------------------------------')
# sp.pprint(t_1in0)
# print('T 2 in 1-------------------------------------------')
# sp.pprint(t_2in1)
# print('T 3 in 2-------------------------------------------')
# sp.pprint(t_3in2)
# print('T 4 in 3-------------------------------------------')
# sp.pprint(t_4in3)
# print('T 4 in 0-------------------------------------------')
# sp.pprint(t_4in0)

omega_0in0 = sp.Matrix([0, 0, 0])
omega_1in1 = calculate_zdh_rotation_speed(omega_0in0, t_1in0, q1, t)
omega_2in2 = calculate_zdh_rotation_speed(omega_1in1, t_2in1, q2, t)
omega_3in3 = calculate_zdh_rotation_speed(omega_2in2, t_3in2, q3, t)
omega_4in4 = calculate_zdh_rotation_speed(omega_3in3, t_4in3)

# print('omega 0 in 0-------------------------------------------')
# sp.pprint(omega_0in0)
# print('omega 1 in 1-------------------------------------------')
# sp.pprint(omega_1in1)
# print('omega 2 in 2-------------------------------------------')
# sp.pprint(omega_2in2)
# print('omega 3 in 3-------------------------------------------')
# sp.pprint(omega_3in3)
# print('omega 4 in 4-------------------------------------------')
# sp.pprint(omega_4in4)

v_0in0 = sp.Matrix([0, 0, 0])
v_1in1 = calculate_zdh_linear_speed(v_0in0, omega_0in0, t_1in0)
v_2in2 = calculate_zdh_linear_speed(v_1in1, omega_1in1, t_2in1)
v_3in3 = calculate_zdh_linear_speed(v_2in2, omega_2in2, t_3in2)
v_4in4 = calculate_zdh_linear_speed(v_3in3, omega_3in3, t_4in3)

# print('v 0 in 0-------------------------------------------')
# sp.pprint(v_0in0)
# print('v 1 in 1-------------------------------------------')
# sp.pprint(v_1in1)
# print('v 2 in 2-------------------------------------------')
# sp.pprint(v_2in2)
# print('v 3 in 3-------------------------------------------')
# sp.pprint(v_3in3)
# print('v 4 in 4-------------------------------------------')
# sp.pprint(v_4in4)

m1 = sp.Symbol('m1')
m2 = sp.Symbol('m2')
m3 = sp.Symbol('m3')
m4 = sp.Symbol('m4')

I1xx = sp.Symbol('I1xx')
I1yy = sp.Symbol('I1yy')
I1zz = sp.Symbol('I1zz')

I2xx = sp.Symbol('I2xx')
I2yy = sp.Symbol('I2yy')
I2zz = sp.Symbol('I2zz')

I3xx = sp.Symbol('I3xx')
I3yy = sp.Symbol('I3yy')
I3zz = sp.Symbol('I3zz')

I4xx = sp.Symbol('I4xx')
I4yy = sp.Symbol('I4yy')
I4zz = sp.Symbol('I4zz')

I1 = sp.Matrix([[I1xx, 0, 0],
                [0, I1yy, 0],
                [0, 0, I1zz]])
I2 = sp.Matrix([[I2xx, 0, 0],
                [0, I2yy, 0],
                [0, 0, I2zz]])
I3 = sp.Matrix([[I3xx, 0, 0],
                [0, I3yy, 0],
                [0, 0, I3zz]])
# I4 = sp.Matrix([[I4xx, 0, 0],
#                 [0, I4yy, 0],
#                 [0, 0, I4zz]])

c1z = sp.Symbol('c1z')
c2x = sp.Symbol('c2x')
c3x = sp.Symbol('c3x')

c1 = sp.Matrix([0, 0, c1z])
c2 = sp.Matrix([c2x, 0, 0])
c3 = sp.Matrix([c3x, 0, 0])
# c4 = sp.Matrix([0, 0, 0])

g = sp.Symbol('g')
gravity = sp.Matrix([0, 0, g])

Ek1 = calculate_kinetic_energy(omega_1in1, v_1in1, I1, c1, m1)
Ek2 = calculate_kinetic_energy(omega_2in2, v_2in2, I2, c2, m2)
Ek3 = calculate_kinetic_energy(omega_3in3, v_3in3, I3, c3, m3)
# Ek4 = calculate_kinetic_energy(omega_4in4, v_4in4, I4, c4, m4)

Ep1 = calculate_potential_energy(m1, gravity, t_1in0, c1)
Ep2 = calculate_potential_energy(m2, gravity, t_2in0, c2)
Ep3 = calculate_potential_energy(m3, gravity, t_3in0, c3)
# Ep4 = calculate_potential_energy(m4, gravity, t_4in0, c4)

# print('kinetic energy 1-------------------------------------------')
# sp.pprint(Ek1)
# print('kinetic energy 2-------------------------------------------')
# sp.pprint(Ek2)
# print('kinetic energy 3-------------------------------------------')
# sp.pprint(Ek3)
# print('kinetic energy 4-------------------------------------------')
# sp.pprint(Ek4)

# print('potential energy 1-------------------------------------------')
# sp.pprint(Ep1)
# print('potential energy 2-------------------------------------------')
# sp.pprint(Ep2)
# print('potential energy 3-------------------------------------------')
# sp.pprint(Ep3)
# print('potential energy 4-------------------------------------------')
# sp.pprint(Ep4)

Ekc = Ek1 + Ek2 + Ek3 # + Ek4
Epc = Ep1 + Ep2 + Ep3 # + Ep4

L = Ekc - Epc

# print('Lagrange ---------------------------------------------------')
# sp.pprint(L)

q1_prim = q1.diff(t)
q2_prim = q2.diff(t)
q3_prim = q3.diff(t)
cords_prim = sp.Matrix([q1_prim, q2_prim, q3_prim])

M = calculate_mass_matrix(L, cords, cords_prim, 3, t)

print('----------chwila prawdy----------')
sp.pprint(M)


C = calculate_coriolis_matrix(M, cords, cords_prim)

# file = open('coriolis_matrix.tex', 'w+')
# file.write(sp.latex(C))
# file.close()


G = calculate_gravitation_vector(Epc, cords, 3)
# sp.pprint(G)

x1 = cords
x2 = x1_prim = cords_prim

# x1' = x2
# M(x1)x2' + C(x1,x2){x2} + G(x1) = u

# M_inversed = M.inv()
# sp.pprint(M_inversed)

# podstawienie zmiennych do modelu
