import math
import sympy as sp


def calculate_kinetic_energy(rotation_speed: sp.Matrix, linear_speed: sp.Matrix,
                             inertial_tensor: sp.Matrix, mass_cords: sp.Matrix, mass):
    valid_shape = rotation_speed.shape == tuple([3, 1]) and linear_speed.shape == tuple([3, 1]) and \
                  inertial_tensor.shape == tuple([3, 3]) and mass_cords.shape == tuple([3, 1])
    if not valid_shape:
        raise ValueError('Not valid matrices shapes')

    kinetic_energy = 0.5 * ((rotation_speed.transpose() * inertial_tensor * rotation_speed) +
                            (mass * linear_speed.transpose() * linear_speed) +
                            (2 * linear_speed.transpose() * (rotation_speed.cross(mass * mass_cords))))
    return kinetic_energy


def calculate_potential_energy(mass, gravity: sp.Matrix, zdh_matrix: sp.Matrix, mass_cords: sp.Matrix):
    valid_shape = gravity.shape == tuple([3, 1]) and zdh_matrix.shape == tuple([4, 4]) and \
                  mass_cords.shape == tuple([3, 1])
    if not valid_shape:
        raise ValueError('Not valid parameters')

    r = zdh_matrix[0:3, 0:3]
    l = zdh_matrix[0:3, -1]
    # sp.pprint(l)

    potential_energy = -mass * gravity.transpose() * (l + r * mass_cords)

    return potential_energy


def calculate_coriolis_matrix(mass_matrix: sp.Matrix, state_coordinates: sp.Matrix, state_coordinates_prims: sp.Matrix):
    coriolis_matrix = sp.ImmutableDenseMatrix(sp.zeros(3))

    for k in range(0, mass_matrix.cols):
        for j in range(0, mass_matrix.rows):
            for i in range(0, state_coordinates.rows):
                print(k, j, i)
                coriolis_matrix[k, j] += (mass_matrix[k, j].diff(state_coordinates[i]) +
                                          mass_matrix[k, i].diff(state_coordinates[j]) -
                                          mass_matrix[i, j].diff(state_coordinates[k])) * state_coordinates_prims[i]
                # coriolis_matrix[k, j] = temp

    return coriolis_matrix


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
        rotation_speed += state_function.diff(t) * sp.Matrix([0, 0, 1])

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
        linear_speed += state_function.diff(t) * sp.Matrix([0, 0, 1])

    return linear_speed


def calculate_zdh_matrix(theta_i, alpha_i_minus1, a_i_minus1, d_i):
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


def s(value):
    if type(type(value)) == sp.function.UndefinedFunction:
        return sp.sin(value)
    elif type(value) == float:
        temp_trigonometry_value = sp.sin(value).evalf()
        if 0.9999 < temp_trigonometry_value:
            return 1.0
        elif -0.0001 < temp_trigonometry_value < 0.0001:
            return 0.0
        elif temp_trigonometry_value < -0.9999:
            return -1.0
        else:
            return temp_trigonometry_value
    else:
        raise ValueError('bad value was passed to the function')


def c(value):
    if type(type(value)) == sp.function.UndefinedFunction:
        return sp.cos(value)
    elif type(value) == float:
        temp_trigonometry_value = sp.cos(value).evalf()
        if 0.9999 < temp_trigonometry_value:
            return 1.0
        elif -0.0001 < temp_trigonometry_value < 0.0001:
            return 0.0
        elif temp_trigonometry_value < -0.9999:
            return -1.0
        else:
            return temp_trigonometry_value
    else:
        return None


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
I4 = sp.Matrix([[I4xx, 0, 0],
                [0, I4yy, 0],
                [0, 0, I4zz]])

c1z = sp.Symbol('c1z')
c2x = sp.Symbol('c2x')
c3x = sp.Symbol('c3x')

c1 = sp.Matrix([0, 0, c1z])
c2 = sp.Matrix([c2x, 0, 0])
c3 = sp.Matrix([c3x, 0, 0])
c4 = sp.Matrix([0, 0, 0])

g = sp.Symbol('g')
gravity = sp.Matrix([0, 0, g])

Ek1 = calculate_kinetic_energy(omega_1in1, v_1in1, I1, c1, m1)
Ek2 = calculate_kinetic_energy(omega_2in2, v_2in2, I2, c2, m2)
Ek3 = calculate_kinetic_energy(omega_3in3, v_3in3, I3, c3, m3)
Ek4 = calculate_kinetic_energy(omega_4in4, v_4in4, I4, c4, m4)

Ep1 = calculate_potential_energy(m1, gravity, t_1in0, c1)
Ep2 = calculate_potential_energy(m2, gravity, t_2in0, c2)
Ep3 = calculate_potential_energy(m3, gravity, t_3in0, c3)
Ep4 = calculate_potential_energy(m4, gravity, t_4in0, c4)

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

Ekc = Ek1 + Ek2 + Ek3 + Ek4
Epc = Ep1 + Ep2 + Ep3 + Ep4

L = Ekc - Epc

# print('Lagrange ---------------------------------------------------')
# sp.pprint(L)

q1_prim = q1.diff(t)
q2_prim = q2.diff(t)
q3_prim = q3.diff(t)
cords_prim = sp.Matrix([q1_prim, q2_prim, q3_prim])

q1_bis = q1_prim.diff(t)
q2_bis = q2_prim.diff(t)
q3_bis = q3_prim.diff(t)

L_q1_prim_t = L.diff(q1_prim).diff(t)
L_q2_prim_t = L.diff(q2_prim).diff(t)
L_q3_prim_t = L.diff(q3_prim).diff(t)

# print('L _q3_prim po t - -----------------')
# sp.pprint(L_q3_prim_t)

L_q1 = L.diff(q1)
L_q2 = L.diff(q2)
L_q3 = L.diff(q3)

# M MATRIX FIRST ROW
m11 = L_q1_prim_t.subs(q2_bis, 0)
m11 = m11.subs(q3_bis, 0)
m11 = m11.subs(q2_prim, 0)
m11 = m11.subs(q3_prim, 0)
m11 = m11.subs(q1_bis, 1)
m11 = m11.subs(q1_prim, 0)

m12 = L_q1_prim_t.subs(q1_bis, 0)
m12 = m12.subs(q3_bis, 0)
m12 = m12.subs(q1_prim, 0)
m12 = m12.subs(q3_prim, 0)
m12 = m12.subs(q2_bis, 1)
m12 = m12.subs(q2_prim, 0)

m13 = L_q1_prim_t.subs(q1_bis, 0)
m13 = m13.subs(q2_bis, 0)
m13 = m13.subs(q1_prim, 0)
m13 = m13.subs(q2_prim, 0)
m13 = m13.subs(q3_bis, 1)
m13 = m13.subs(q3_prim, 0)

# M MATRIX SECOND ROW
m21 = L_q2_prim_t.subs(q2_bis, 0)
m21 = m21.subs(q3_bis, 0)
m21 = m21.subs(q2_prim, 0)
m21 = m21.subs(q3_prim, 0)
m21 = m21.subs(q1_bis, 1)
m21 = m21.subs(q2_prim, 0)

m22 = L_q2_prim_t.subs(q1_bis, 0)
m22 = m22.subs(q3_bis, 0)
m22 = m22.subs(q1_prim, 0)
m22 = m22.subs(q3_prim, 0)
m22 = m22.subs(q2_bis, 1)
m22 = m22.subs(q2_prim, 0)

m23 = L_q2_prim_t.subs(q2_bis, 0)
m23 = m23.subs(q3_bis, 0)
m23 = m23.subs(q2_prim, 0)
m23 = m23.subs(q3_prim, 0)
m23 = m23.subs(q1_bis, 1)
m23 = m23.subs(q1_prim, 0)

# M Matrix THIRD ROW
m31 = L_q3_prim_t.subs(q2_bis, 0)
m31 = m31.subs(q3_bis, 0)
m31 = m31.subs(q2_prim, 0)
m31 = m31.subs(q3_prim, 0)
m31 = m31.subs(q1_bis, 1)
m31 = m31.subs(q1_prim, 0)
# print('----------m31----------')
# sp.pprint(m31)

m32 = L_q3_prim_t.subs(q1_bis, 0)
m32 = m32.subs(q3_bis, 0)
m32 = m32.subs(q1_prim, 0)
m32 = m32.subs(q3_prim, 0)
m32 = m32.subs(q2_bis, 1)
m32 = m32.subs(q2_prim, 0)

m33 = L_q3_prim_t.subs(q1_bis, 0)
m33 = m33.subs(q2_bis, 0)
m33 = m33.subs(q1_prim, 0)
m33 = m33.subs(q2_prim, 0)
m33 = m33.subs(q3_bis, 1)
m33 = m33.subs(q3_prim, 0)
# print('----------m33----------')
# sp.pprint(m33)

M = sp.Matrix([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])
# print('----------chwila prawdy----------')
# sp.pprint(M)


sp.pprint(M[0, 0].diff(cords[0]))
# sp.pprint(cords)
calculate_coriolis_matrix(M, cords, cords_prim)
