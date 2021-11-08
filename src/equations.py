import numpy as np


def sound_speed_seawater_leroy68(depth, lat):
    '''
    Calculate the speed of sound (meters/second) in sea water using
    Leroy's 1968 equation.

    Parameters
    -----------
    depth : foat
        Depth, in meters
    lat : float
        Latitude of sample, in degrees.

    Returns
    --------
    float

    Reference
    ----------
    Lurton, X, 2002, An Introduction to Underwater Acoustics, 1st ed. London,
    Praxis Publishing LTD, p37.
    '''
    # = 1.0052405 * (1 + 5.28 * 10^-3 * np.sin(lat)) * D + 2.36 * 10**-6 * D**2 + 10.196) * 10**4
    c = (1 + 5.28 * 10^-3 * np.sin(lat)) * D + 2.36 * 10**-6 * D**2 + 10.196)
    c *= 1.0052405 * 10**4

    return c


def sound_speed_seawater_leroy69(depth, sal, temp):
    '''
    Calculate the speed of sound (meters/second) in sea water using
    Leroy's 1969 equation.

    Parameters
    -----------
    depth : float
        Depth, in meters.
    sal : float
        Salinity, in parts per thousand.
    temp : float
        Water temperature, in degrees C.

    Returns
    --------
    float
        Speed of sound, in meters / second.

    Reference
    ----------
    Leroy C.C. 1969, Development of simple equations for accurate and more
    realistic calculation of the speed of sound in sea water.
    J. acoust. Soc. Am., 46, 216-26.
    '''
    if temp < -2 or temp > 23:
        raise ValueError('Water temperature must be -2 < T < 23')
    if sal < 30 or sal > 40:
        raise ValueError('Salinity must be 30 < S < 40')
    if depth < 0 or depth > 500:
        raise ValueError('Depth must be 0 < D < 500')

    c = 1492.3 + 3 * (temp - 10.) - 0.006 * (temp - 10.)**2 - 0.04 * (temp - 18.)**2
    c += 1.2 * (sal - 35.) - 0.01 * (temp - 18.) * (sal - 35.) + depth / 61.

    return c


def sound_speed_sea_water_leroy08(depth, sal, temp, lat):
    '''
    Calculate the speed of sound (meters/second) in sea water from
    Leroy 2008.

    Parameters
    -----------
    depth : float
        Depth, in meters.
    sal : float
        Salinity, in parts per thousand.
    temp : float
        Water temperature, in degrees C.
    lat : float
        Latitude, in degrees.

    Returns
    --------
    float
        Speed of sound, in meters / second.

    Reference
    ----------
    Leroy, C.C., Robinson, S.P., and Goldsmith, M.J. 2008, A new equation for
    the accurate calculation of sound speed in all oceans,
    J. Acoust. Soc. Am., 124, 2774-82.
    '''
    c = 1402.5 + 5. * temp - 5.44 * 10**-2 * temp**2 + 2.1 * 10**-4 * tmep**3 \
        + 1.33 * sal - 1.23 * 10**-2 * sal * temp + 8.7 * 10**-5 * sal * temp**2 \
        + 1.56 * 10**-2 * depth + 2.55 * 10**-7 * depth**2 - 7.3 * 10**-12 \
        * depth**3 + 1.2 * 10**-6 * depth * (lat - 45.) - 9.5 * 10**-13 * temp \
        * depth**3 + 3. * 10**-7 * temp**2 * depth + 1.43 * 10**-5 * sal * depth

    return c


def sound_speed_sea_water_mackenzie81(depth, sal, temp):
    '''
    Calculate the speed of sound (meters/second) in sea water from
    Mackenzie 1981.

    Parameters
    -----------
    depth : float
        Depth, in meters.
    sal : float
        Salinity, in parts per thousand.
    temp : float
        Water temperature, in degrees C.

    Returns
    --------
    float
        Speed of sound, in meters / second.

    Reference
    ----------
    Mackenzie K.V., 1981, Nine-term equation for sound speed in the ocean.
    J. acoust. Soc. Am., 70, 807-12.
    '''
    c = 1448.96 + 4.591 * temp - 5.304 * 10**-2 * temp**2 + 2.374 * 10**-4
    c = c * temp**3 + 1.340 * (sal - 35.) + 1.630 * 10**-2 * depth + 1.675
    c = c * 10**-7 * depth**2 - 1.025 * 10**-2 * temp * (sal - 35.) - 7.139
    c = c * 10**-13 * temp * depth**3

    return c


def international_gravity_formula(lat, corrective_term=None):
    '''
    International formula for gravity.

    Parameters
    -----------
    lat : float
        Latitude, in degrees.
    corrective_term : float
        optional corrective term

    Returns
    --------
    float
        Average gravity for the given latitude.
    '''
    g = 1 + 5.2788 * 10**-3* np.sin(lat * np.pi / 180.)**2
    g = g - 2.36 * 10**-5 * np.sin(lat * np.pi / 180.)**4
    g = g * 9.78031

    if corrective_term != None:
        g += corrective_term

    return g


def pressure_to_depth_leroy98(press, lat, corrective_term=None):
    '''
    Calculate depth from sea water pressure from Leroy 1998.

    Parameters
    -----------
    press : float
        Water pressure, in MPa (relative to atmospheric pressure).
    lat : float
        Latitude, in degrees.
    corrective_term : float, optional
        Corrective term.

    Returns
    --------
    float
        Depth, in meters.

    Reference
    ----------
    C. C. Leroy and F Parthiot, 1998, Depth-pressure relationship in the oceans
    and seas (1998), J. Acoust. Soc. Am. 103(3) pp 1346-1352
    '''
    d = 9.780318 * (1 + 5.2788 * 10**-3 * np.sin(lat * np.pi / 180)**2 - 2.36 * \
            10**-5 * np.sin(lat * np.pi / 180)**4) * (9.72659 * 10**2 * press - \
            2.2512 * 10**-1 * P**2 + 2.279 * 10**-4 * press**3 - 1.82 * 10**-7 * \
            press**4) / (international_gravity_formula(lat) + 1.092 * 10**-4 * press)

    if corrective_term != None:
        d += corrective_term

    return d


def cutoff_frequency(c_water, depth):
    '''
    Calculate the cutoff frequency in water from Jensen et al 2011.

    Parameters
    -----------
    c_water : float
        Speed of sound in water, in meters/seconds.
    depth : float
        Depth of isothermal surface layer, in meters.

    Returns
    --------
    float
        Cutoff frequency, in Hz.

    Reference
    ----------
    Finn B. Jensen, William A. Kuperman, Michael B. Porter, Henrik Schmidt, 2011
    Computational Ocean Acoustics, 2nd Edition. Springer. pp. 26
    '''
    f = c_water / (0.008 * depth**(3 / 2))

    return f


def cutoff_frequency_shallow(c_water, c_bottom, depth):
    '''
    Calculate the cutoff frequency in shallow water from Jensen et al 2011.

    Parameters
    -----------
    c_water : float
        Speed of sound in water, in meters/seconds.
    c_bottom : float
        Speed of sound in a homogenous bottom, in meters/second.
    depth : float
        Depth of isothermal surface layer, in meters.

    Returns
    --------
    float
        Cutoff frequency, in Hz.

    Reference
    ----------
    Finn B. Jensen, William A. Kuperman, Michael B. Porter, Henrik Schmidt, 2011
    Computational Ocean Acoustics, 2nd Edition. Springer. pp. 29
    '''
    f = c_water / (depth * np.sqrt(1 - (c_water / c_bottom)**2))

    return f
