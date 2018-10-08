import math
from robotlib import *


class Robot:
    g = sp.Symbol('g')
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
