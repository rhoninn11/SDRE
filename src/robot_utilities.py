import sympy as sp


def s(value):
    if type(type(value)) == sp.function.UndefinedFunction or type(value) == sp.Symbol:
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
    if type(type(value)) == sp.function.UndefinedFunction or type(value) == sp.Symbol:
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