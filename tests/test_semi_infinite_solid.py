"""The semi-infinite solid, what a wonderful specimen

###############################################################################
# email:    noahgill409@gmail.com                                             #
# package:  cheme                                                             #
###############################################################################

Geankoplis example 5.3-1

The depth in the soil of the earth at which freezing temperatures penetrate is
often of importance in agriculture and construction. During a certain fall
day, the temperature in the earth is a constant at 15.6 degC to a depth of
several meters. A cold wave suddenly reduces the air temperature from 15.6 to
-17.8. The convective coefficient above the soil is 11.36 W/m^2*K. The soil
properties can be assumed as alpha = 4.65e-7 m^2/s and k = 0.865 W/m*K.
Neglect any latent heat effects.
"""

import math

from cheme.transfer.semi_infinite_solid import semi_infinite_solid


def test_semi_infinite_solid():
    T0 = 15.6 + 273.15      # K
    T1 = -17.8 + 273.15     # K
    h = 11.36               # W / m^2 * K
    alpha = 4.65E-7         # m^2 / s
    k = 0.865               # W / m * K

    # part a
    time = 5 * 3600         # seconds
    x = 0                   # m (surface)

    T = semi_infinite_solid(time, T0, T1, h, x, alpha, k)

    assert math.fabs(T - 267.982) < 0.001
