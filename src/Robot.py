import math
from robotlib import *
from sympy.printing.latex import print_latex


class Robot:
    # g = sp.Symbol('g')
    g = -10
    t = 0

    gravity = sp.Matrix([0, 0, g])

    links_number = 0
    cords_list = []
    cords_prim_list = []

    links = []

    transformation_matrices = []
    transformation_matrices_in_base = []

    omegas = []
    velocities = []

    Ekc = sp.Matrix([0])
    Epc = sp.Matrix([0])
    L = 0

    M = 0
    C = 0
    G = 0

    cords = 0
    cords_prim = 0

    def __init__(self, t: sp.Symbol):
        self.t = t

    def add_link(self, link: Link):
        if link.type == LinkType.linear:
            self.cords_list.append(link.d)
        elif link.type == LinkType.rotational:
            self.cords_list.append(link.theta)

        self.links.append(link)
        if not link.type == LinkType.none:
            self.links_number += 1

    def calculate_dynamics(self):
        self.define_cords()
        self.calculate_transformations_matrices()
        self.calculate_tranformation_matrices_in_base()
        self.calculate_velocities()
        self.calculate_lagrangian()

        self.M = calculate_mass_matrix(self.L, self.cords, self.cords_prim, self.links_number, self.t)
        self.C = calculate_coriolis_matrix(self.M, self.cords, self.cords_prim, self.links_number)
        self.G = calculate_gravitation_vector(self.Epc, self.cords, self.links_number)

    def define_cords(self):
        for i in range(0, self.links_number):
            if not self.links[i].type == LinkType.none:
                self.cords_prim_list.append(self.cords_list[i].diff(self.t))

        self.cords = sp.Matrix(self.cords_list)
        self.cords_prim = sp.Matrix(self.cords_prim_list)

    def calculate_transformations_matrices(self):
        #  macierze T (T[0] = T1in0, T[1] = T_2in1)------------------
        self.transformation_matrices = []
        for link in self.links:
            self.transformation_matrices.append(calculate_zdh_matrix(link))

    def calculate_tranformation_matrices_in_base(self):
        #  macierze w układzie 0--------------------------
        self.transformation_matrices_in_base = []
        for i in range(0, len(self.transformation_matrices)):
            temp = 1
            for j in range(0, i + 1):
                temp *= self.transformation_matrices[j]
            self.transformation_matrices_in_base.append(temp)

    def calculate_velocities(self):
        self.calculate_rotation_velocities()
        self.calculate_linear_velocities()

    def calculate_rotation_velocities(self):
        #  prędkości kątowe-----------------------------
        self.omegas = []  # omegas[0] = omega_0in0 ...
        self.omegas.append(sp.Matrix([0, 0, 0]))
        for i in range(0, len(self.transformation_matrices)):
            if self.links[i].type == LinkType.rotational:
                self.omegas.append(calculate_zdh_rotation_speed(self.omegas[i], self.transformation_matrices[i],
                                                                self.cords_list[i], self.t))
            else:
                self.omegas.append(calculate_zdh_rotation_speed(self.omegas[i], self.transformation_matrices[i]))

    def calculate_linear_velocities(self):
        #  prędkości liniowe---------------------------
        self.velocities = []
        self.velocities.append(sp.Matrix([0, 0, 0]))
        for i in range(0, len(self.transformation_matrices)):
            if self.links[i].type == LinkType.linear:
                self.velocities.append(calculate_zdh_linear_speed(self.velocities[i], self.omegas[i],
                                                                  self.transformation_matrices[i],
                                                                  self.cords_list[i], self.t))
            else:
                self.velocities.append(calculate_zdh_linear_speed(self.velocities[i], self.omegas[i],
                                                                  self.transformation_matrices[i]))

    def calculate_lagrangian(self):
        kinetic_energy = []
        potential_energy = []
        for i in range(0, len(self.cords_list)):
            kinetic_energy.append(calculate_kinetic_energy(self.omegas[i + 1], self.velocities[i + 1], self.links[i]))
            potential_energy.append(
                calculate_potential_energy(self.gravity, self.transformation_matrices_in_base[i], self.links[i]))
            self.Ekc += kinetic_energy[i]
            self.Epc += potential_energy[i]
        self.L = self.Ekc - self.Epc

    def print_dynamics(self):
        print('-------------------M - mass matrix------------------------')
        sp.pprint(self.M)
        print('-------------------C - corioliss matrix-------------------')
        sp.pprint(self.C)
        print('-------------------G - gravitation vector-----------------')
        sp.pprint(self.G)

    def print_transofrmation_matrices(self):
        for i in range(0, len(self.transformation_matrices)):
            print('-------------------T matrix ' + str(i + 1) + ' in ' + str(i) + '------------------------')
            sp.pprint(self.transformation_matrices[i])

    def print_transofrmation_matrices_in_base(self):
        for i in range(0, len(self.transformation_matrices_in_base)):
            print('-------------------T matrix ' + str(i + 1) + ' in ' + str(0) + '------------------------')
            sp.pprint(self.transformation_matrices_in_base[i])

    def print_veocitioes(self):
        for i in range(0, len(self.omegas)):
            print('------------------ omega ' + str(i) + ' in ' + str(i) + '------------------------')
            sp.pprint(self.omegas[i])
        for i in range(0, len(self.velocities)):
            print('-------------------v ' + str(i) + ' in ' + str(i) + '------------------------')
            sp.pprint(self.velocities[i])

    def print_to_latex(self):
        f = open('robot.tex', '+w')
        f.write('\\documentclass[10]{article}\n')
        f.write('\\usepackage{amsmath}\n')
        f.write('\\usepackage[paperheight = 11cm, paperwidth = 500cm]{geometry}\n')
        f.write('\\begin{document}\n')
        f.write('------------------------Mass matrix---------------------------\\\\\n')
        f.write('\\\\$')
        f.write(sp.latex(self.M))
        f.write('$\\\\\n')
        f.write('------------------------Corioliss matrix---------------------------\\\\\n')
        f.write('\\\\$')
        f.write(sp.latex(self.C))
        f.write('$\\\\\n')
        f.write('------------------------Gravitation vector---------------------------\\\\\n')
        f.write('\\\\$')
        f.write(sp.latex(self.G))
        f.write('$\\\\\n')
        f.write('\\end{document}\n')
        f.close()

    def print_M_elements_to_latex(self):
        for i in range(0, self.links_number):
            for j in range(0, self.links_number):
                print('M-element', i + 1, j + 1, '--------------------------')
                sting_latex = sp.latex(self.M[i, j])
                print(type(sting_latex))
                sting_latex.replace('\left', '')
                sting_latex.replace('\right', '')
                print(sting_latex)
                print('------------------------------------------')

    def print_C_elements_to_latex(self):
        for i in range(0, self.links_number):
            for j in range(0, self.links_number):
                print('C-element', i + 1, j + 1, '--------------------------')
                sting_latex = sp.latex(self.C[i, j])
                sting_latex.replace('\left', '')
                sting_latex.replace('\right', '')
                print(sting_latex)
                print('------------------------------------------')

    def print_G_elements_to_latex(self):
        for i in range(0, self.links_number):
            print('G-element', i + 1, '--------------------------')
            sting_latex = sp.latex(self.G[i])
            sting_latex.replace('\left', '')
            sting_latex.replace('\right', '')
            print(sting_latex)
            print('----------------------------------------')

    def insert_cords_values(self):

        x1 = [sp.Function('x11')(self.t),sp.Function('x12')(self.t),sp.Function('x13')(self.t)]
        x2 =  x1_prim = [sp.Function('x21')(self.t),sp.Function('x22')(self.t),sp.Function('x23')(self.t)]
        M_obliczeniowa = 0
        C_obliczeniowa = 0
        G_obliczeniowa = 0
        for i in range(len(self.cords_list)):
            if i is 0:
                M_obliczeniowa = self.M.subs(self.cords_list[i],x1[i])
                C_obliczeniowa = self.C.subs(self.cords_prim_list[i],x2[i])
                C_obliczeniowa = C_obliczeniowa.subs(self.cords_list[i],x1[i])
                G_obliczeniowa = self.G.subs(self.cords_list[i],x1[i])
            else:
                M_obliczeniowa = M_obliczeniowa.subs(self.cords_list[i], x1[i])
                C_obliczeniowa = C_obliczeniowa.subs(self.cords_prim_list[i],x2[i])
                C_obliczeniowa = C_obliczeniowa.subs(self.cords_list[i],x1[i])
                G_obliczeniowa = G_obliczeniowa.subs(self.cords_list[i], x1[i])

        #to mamy macxierz mas i C i G gotową

        sp.pprint(M_obliczeniowa)
        print('------')
        sp.pprint(C_obliczeniowa)
        print('------')
        sp.pprint(G_obliczeniowa)
        print('------')

        x1_values = [0.5,0.6,0.8]
        x2_values = [1.2,0.1,0.6]

        for i in range(len(self.cords_list)):
            M_obliczeniowa = M_obliczeniowa.subs(x1[i],x1_values[i])
            C_obliczeniowa = C_obliczeniowa.subs(x2[i],x2_values[i])
            C_obliczeniowa = C_obliczeniowa.subs(x1[i],x1_values[i])
            G_obliczeniowa = G_obliczeniowa.subs(x1[i],x1_values[i])

        sp.pprint(M_obliczeniowa)
        print('------')
        sp.pprint(C_obliczeniowa)
        print('------')
        sp.pprint(G_obliczeniowa)
        print('------')
