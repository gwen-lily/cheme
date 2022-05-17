"""The semi-infinite solid, what a wonderful specimen

###############################################################################
# email:    noahgill409@gmail.com                                             #
# package:  cheme                                                             #
###############################################################################
"""


import math

import numpy as np
from scipy.optimize import fsolve


def wumbo(
    T: float,
    t: float,
    T0: float,
    T1: float,
    h: float | int,
    x: float,
    alpha: float,
    k: float | int
) -> float:
    """unsteady-state conduction equation in a semi-infinite solid

    See Geankoplis 5.3B for a full explanation. Originally, the temperature in
    the solid is uniform at T0. At time t = 0, the solid is suddenly exposed
    to or immersed in a large mass of ambient fluid at tmemperature T1, which
    is constant. The convection coefficient h in W/m^2*K or btu/h*ft^2*degF is
    present and is constant; i.e., a surface resistance is present. Hence, the
    temperature TS at the surface is not the same as T1.


    Parameters
    ----------
    T : float
        The temperature at the relevant position and time
    time : float
        time in seconds
    T0 : float
        The original uniform temperature of the solid in Kelvin
    T1 : float
        The temperature of the ambient fluid in Kelvin
    h : float
        _description_
    x : float
        the distance into the solid from the surface in SI units in meters
    alpha : float
        alpha = k/rho*c_p in m^2/s.
    k : float
        _description_

    Returns
    -------
    float
        temperature
    """

    fraction_of_change = (T - T0) / (T1 - T0)

    # shorthand parameters
    m = math.sqrt(alpha*t)
    x_2m = x / (2 * m)
    hm_k = h * m / k
    exp_term = hm_k * (x / m + hm_k)

    zero = -1 * fraction_of_change + math.erfc(x_2m) - \
        math.exp(exp_term) * math.erfc(x_2m + hm_k)

    return zero


def solve_for_T(*args, **kwargs) -> float:
    return wumbo(*args, **kwargs)


def solve_for_t(
    t: float,
    T: float,
    T0: float,
    T1: float,
    h: float | int,
    x: float,
    alpha: float,
    k: float | int
) -> float:
    return wumbo(T, t, T0, T1, h, x, alpha, k)


def find_temperature(
    t: float,
    T0: float,
    T1: float,
    h: float,
    x: float,
    alpha: float,
    k: float,
    initial_guess: float | None = None,
) -> float:
    """unsteady-state conduction equation in a semi-infinite solid

    See Geankoplis 5.3B for a full explanation.


    Parameters
    ----------
    t : float
        The time in seconds
    T0 : float
        The original uniform temperature of the solid in Kelvin
    T1 : float
        The temperature of the ambient fluid in Kelvin
    h : float
        convective coefficient in W/m^2*K
    x : float
        the distance into the solid from the surface in SI units in meters
    alpha : float
        solid property: alpha = k/rho*c_p in m^2/s.
    k : float
        solid property

    Returns
    -------
    float
        time in seconds
    """

    if initial_guess is None:
        initial_guess = (T0 + T1) / 2

    soln = fsolve(solve_for_T, initial_guess, (t, T0, T1, h, x, alpha, k))

    assert isinstance(soln, np.ndarray)
    return soln[0]


def find_time(
    T: float,
    T0: float,
    T1: float,
    h: float,
    x: float,
    alpha: float,
    k: float,
    initial_guess: float | None = None
) -> float:
    """unsteady-state conduction equation in a semi-infinite solid

    See Geankoplis 5.3B for a full explanation.


    Parameters
    ----------
    T : float
        The temperature at the specified point in Kelvin
    T0 : float
        The original uniform temperature of the solid in Kelvin
    T1 : float
        The temperature of the ambient fluid in Kelvin
    h : float
        convective coefficient in W/m^2*K
    x : float
        the distance into the solid from the surface in SI units in meters
    alpha : float
        solid property: alpha = k/rho*c_p in m^2/s.
    k : float
        solid property

    Returns
    -------
    float
        time in seconds
    """

    if initial_guess is None:
        initial_guess = (k / h) ** 2 / alpha

    soln = fsolve(solve_for_t, initial_guess, (T, T0, T1, h, x, alpha, k))

    assert isinstance(soln, np.ndarray)
    return soln[0]
