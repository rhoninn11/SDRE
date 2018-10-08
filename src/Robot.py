import math
from robotlib import *


class Robot:
    g = sp.Symbol('g')
    gravity = sp.Matrix([0, 0, g])
    def __init__(self):
        self.i = 3

        gravity = sp.Matrix([0, 0, g])

        Ek1 = calculate_kinetic_energy(omega_1in1, v_1in1, I1, c1, m1)
        Ek2 = calculate_kinetic_energy(omega_2in2, v_2in2, I2, c2, m2)
        Ek3 = calculate_kinetic_energy(omega_3in3, v_3in3, I3, c3, m3)
        # Ek4 = calculate_kinetic_energy(omega_4in4, v_4in4, I4, c4, m4)

        Ep1 = calculate_potential_energy(m1, gravity, t_1in0, c1)
        Ep2 = calculate_potential_energy(m2, gravity, t_2in0, c2)
        Ep3 = calculate_potential_energy(m3, gravity, t_3in0, c3)
        # Ep4 = calculate_potential_energy(m4, gravity, t_4in0, c4)

        Ekc = Ek1 + Ek2 + Ek3  # + Ek4
        Epc = Ep1 + Ep2 + Ep3  # + Ep4

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


    def add_link(self):


