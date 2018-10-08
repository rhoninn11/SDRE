import math
from robotlib import *
import sympy as sp
import sympy.printing.latex as spltx

t = sp.Symbol('t')
links_number = 3
cords_list = []
cords_prim_list = []

for i in range(0, links_number):
    cords_list.append(sp.Function('q'+(i+1)(t)))
    cords_prim_list.append(cords_list[i].diff(t))

#  q1 = sp.Function('q1')(t)
#  q2 = sp.Function('q2')(t)
#  q3 = sp.Function('q3')(t)

cords = sp.Matrix(cords_list)
cords_prim = sp.Matrix(cords_prim_list)

#  definicja przegubów
links = []
links.append(cords[0], 0.0, 0.0, 0.689, 1)
links.append(cords[1], -math.pi/2, 0.0, 0.0, 2)
links.append(cords[2], 0.0, 0.125, 0.0, 3)
links.append(0.0, 0.0, 0.248, 0.0, 4)


#  macierze T (T[0] = T1in0, T[1] = T_2in1)------------------
transformation_matrices = []
for link in links:
    transformation_matrices.append(calculate_zdh_matrix(link))

#  macierze w układzie 0--------------------------
transformation_matrices_in_base = []
for i in range(0, len(transformation_matrices)):
    temp =1
    for j in range(0, i+1):
        temp *= transformation_matrices[j]
    transformation_matrices_in_base.append(temp)

#  prędkości kątowe-----------------------------
omegas = [] # omegas[0] = omega_0in0 ...
omegas.append(sp.Matrix([0, 0, 0]))
for i in range(0, len(transformation_matrices)):
    if links[i].type == LinkType.rotational:
        omegas.append(calculate_zdh_rotation_speed(omegas[i], transformation_matrices[i], cords_list[i], t))
    else:
        omegas.append(calculate_zdh_rotation_speed(omegas[i], transformation_matrices[i]))

#  prędkości liniowe---------------------------
velocities = []
velocities.append(sp.Matrix([0, 0, 0]))
for i in range(0, len(transformation_matrices)):
    if links[i].type == LinkType.linear:
        velocities.append(calculate_zdh_linear_speed(velocities[i], omegas[i], transformation_matrices[i],
                                                     cords_list[i], t))
    else:
        velocities.append(calculate_zdh_linear_speed(velocities[i], omegas[i], transformation_matrices[i]))

# c1z = sp.Symbol('c1z')
# c2x = sp.Symbol('c2x')
# c3x = sp.Symbol('c3x')
#
# c1 = sp.Matrix([0, 0, c1z])
# c2 = sp.Matrix([c2x, 0, 0])
# c3 = sp.Matrix([c3x, 0, 0])
# # c4 = sp.Matrix([0, 0, 0])

g = sp.Symbol('g')
gravity = sp.Matrix([0, 0, g])


# energie kinetyczne ---------------------

kinetic_energy = []
potential_energy = []
Ekc = 0
Epc = 0
for i in range(0, len(cords_list)):
    kinetic_energy.append(calculate_kinetic_energy(omegas[i+1], velocities[i+1], links[i]))
    potential_energy.append(calculate_potential_energy(gravity, transformation_matrices_in_base[i], links[i]))
    Ekc += kinetic_energy[i]
    Epc += potential_energy[i]

L = Ekc - Epc

M = calculate_mass_matrix(L, cords, cords_prim, 3, t)


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
