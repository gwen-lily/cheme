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
    time: float,
    T0: float,
    T1: float,
    h: float,
    x: float,
    alpha: float,
    k: float
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
    m = math.sqrt(alpha*time)
    x_2m = x / (2 * m)
    hm_k = h * m / k
    _term0 = hm_k * (x / m + hm_k)

    zero = -1 * fraction_of_change + math.erfc(x_2m) - \
        math.exp(_term0) * math.erfc(x_2m + hm_k)

    return zero


def semi_infinite_solid(
    time: float,
    T0: float,
    T1: float,
    h: float,
    x: float,
    alpha: float,
    k: float
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
    time : float
        time in seconds
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
        temperature
    """

    _T0 = (T0 + T1) / 2
    _T0 = 267.97

    soln = fsolve(wumbo, _T0, (time, T0, T1, h, x, alpha, k))

    assert isinstance(soln, np.ndarray)
    return soln[0]
