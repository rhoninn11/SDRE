import sympy as sp
from robot_utilities import c, s
from Link import *


#  lagrangian function
#  matrices T0in1 T0in2 ...


def calculate_zdh_matrix(link: Link):
    theta_i = link.theta
    alpha_i_minus1 = link.alpha
    d_i = link.d
    a_i_minus1 = link.a

    t_1_1 = c(theta_i)
    t_1_2 = -s(theta_i)
    t_1_3 = 0.0
    t_1_4 = a_i_minus1

    t_2_1 = s(theta_i) * c(alpha_i_minus1)
    t_2_2 = c(theta_i) * c(alpha_i_minus1)
    t_2_3 = -s(alpha_i_minus1)
    t_2_4 = -d_i * s(alpha_i_minus1)

    t_3_1 = s(theta_i) * s(alpha_i_minus1)
    t_3_2 = c(theta_i) * s(alpha_i_minus1)
    t_3_3 = c(alpha_i_minus1)
    t_3_4 = d_i * c(alpha_i_minus1)

    t_4_1 = 0.0
    t_4_2 = 0.0
    t_4_3 = 0.0
    t_4_4 = 1.0

    zdh_matrix = sp.Matrix([[t_1_1, t_1_2, t_1_3, t_1_4],
                            [t_2_1, t_2_2, t_2_3, t_2_4],
                            [t_3_1, t_3_2, t_3_3, t_3_4],
                            [t_4_1, t_4_2, t_4_3, t_4_4]])

    return zdh_matrix


def calculate_zdh_rotation_speed(previous_rotation_speed: sp.Matrix, zdh_matrix: sp.Matrix,
                                 state_function=None, time_symbol=None) -> sp.Matrix:
    valid_shape = previous_rotation_speed.shape == tuple([3, 1]) and zdh_matrix.shape == tuple([4, 4])
    if not valid_shape:
        raise ValueError('Not valid matrices shapes')

    r = zdh_matrix[0:3, 0:3]
    rotation_speed = r.transpose() * previous_rotation_speed

    valid_symbols = type(type(state_function)) is sp.function.UndefinedFunction and \
                    type(time_symbol) is sp.Symbol
    if valid_symbols:
        rotation_speed += state_function.diff(time_symbol) * sp.Matrix([0, 0, 1])

    return rotation_speed


def calculate_zdh_linear_speed(previous_linear_speed: sp.Matrix, previous_rotation_speed: sp.Matrix,
                               zdh_matrix: sp.Matrix,
                               state_function=None, time_symbol=None) -> sp.Matrix:
    valid_shape = previous_rotation_speed.shape == tuple([3, 1]) and previous_linear_speed.shape == tuple(
        [3, 1]) and zdh_matrix.shape == tuple([4, 4])
    if not valid_shape:
        raise ValueError('Not valid matrices shapes')

    r = zdh_matrix[0:3, 0:3]
    l = zdh_matrix[0:3, -1]
    linear_speed = r.transpose() * (previous_linear_speed + previous_rotation_speed.cross(l))

    valid_symbols = type(type(state_function)) is sp.function.UndefinedFunction and \
                    type(time_symbol) is sp.Symbol
    if valid_symbols:
        linear_speed += state_function.diff(time_symbol) * sp.Matrix([0, 0, 1])

    return linear_speed


def calculate_kinetic_energy(rotation_speed: sp.Matrix, linear_speed: sp.Matrix,
                             link: Link):
    mass_cords = link.mass_center
    inertial_tensor = link.inertial_tensor
    mass = link.mass
    valid_shape = rotation_speed.shape == tuple([3, 1]) and linear_speed.shape == tuple([3, 1]) and \
                  inertial_tensor.shape == tuple([3, 3]) and mass_cords.shape == tuple([3, 1])
    if not valid_shape:
        raise ValueError('Not valid matrices shapes')

    kinetic_energy = 0.5 * ((rotation_speed.transpose() * inertial_tensor * rotation_speed) +
                            (mass * linear_speed.transpose() * linear_speed) +
                            (2 * linear_speed.transpose() * (rotation_speed.cross(mass * mass_cords))))
    return kinetic_energy


def calculate_potential_energy(gravity: sp.Matrix, zdh_matrix: sp.Matrix, link: Link):
    mass_cords = link.mass_center
    mass = link.mass
    valid_shape = gravity.shape == tuple([3, 1]) and zdh_matrix.shape == tuple([4, 4]) and \
                  mass_cords.shape == tuple([3, 1])
    if not valid_shape:
        raise ValueError('Not valid parameters')

    r = zdh_matrix[0:3, 0:3]
    l = zdh_matrix[0:3, -1]
    # sp.pprint(l)

    potential_energy = -mass * gravity.transpose() * (l + r * mass_cords)

    return potential_energy


def calculate_mass_matrix(lagrangian, state_coordinates: sp.Matrix, state_coodinates_prims: sp.Matrix,
                          count: int, t):
    state_coordinates_bis = sp.zeros(count, 1)
    lagrangian_diff_by_cords_prims_and_t = sp.zeros(count, 1)
    lagrangian_diff_by_cords = sp.zeros(count, 1)
    mass_matrix = sp.zeros(count)

    valid_shape = state_coordinates.shape == tuple([count, 1]) and state_coodinates_prims.shape == tuple([count, 1])
    if not valid_shape:
        raise ValueError('Not valid parameters')

    for i in range(0, count):
        state_coordinates_bis[i] = state_coodinates_prims[i].diff(t)
        lagrangian_diff_by_cords_prims_and_t[i] = lagrangian.diff(state_coodinates_prims[i]).diff(t)
        lagrangian_diff_by_cords[i] = lagrangian.diff(state_coordinates[i])

    for i in range(0, count):
        for j in range(0, count):
            temp = lagrangian_diff_by_cords_prims_and_t[i].subs(state_coordinates_bis[j], 1)
            temp = temp.subs(state_coodinates_prims[j], 0)
            for k in range(0, count):
                if not k == j:
                    temp = temp.subs(state_coordinates_bis[k], 0)
                    temp = temp.subs(state_coodinates_prims[k], 0)
            mass_matrix[i, j] = temp

    return mass_matrix


def calculate_coriolis_matrix(mass_matrix: sp.Matrix, state_coordinates: sp.Matrix, state_coordinates_prims: sp.Matrix,
                              count: int):
    coriolis_matrix = sp.zeros(count)
    for k in range(0, mass_matrix.cols):
        for j in range(0, mass_matrix.rows):
            for i in range(0, state_coordinates.rows):
                temp_matrix = sp.zeros(count)
                temp_matrix[k, j] = 0.5 * (mass_matrix[k, j].diff(state_coordinates[i]) +
                                           mass_matrix[k, i].diff(state_coordinates[j]) -
                                           mass_matrix[i, j].diff(state_coordinates[k])) * state_coordinates_prims[i]
                coriolis_matrix += temp_matrix
    return coriolis_matrix


def calculate_gravitation_vector(potential_energy, state_coordinates: sp.Matrix, count: int):
    gravitation_vector = sp.zeros(count, 1)

    valid_shape = state_coordinates.shape == tuple([count, 1])
    if not valid_shape:
        raise ValueError('Not valid matrices shapes')

    for i in range(0, count):
        gravitation_vector[i] = potential_energy.diff(state_coordinates[i])

    return gravitation_vector
