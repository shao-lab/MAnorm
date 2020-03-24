"""
manorm.stats
------------

This module contains mathmatical and statistical functions used in MAnorm.
"""

from math import exp, log


def xy_to_ma(x, y):
    r"""Calculate (M, A) value with given read counts/densities of two samples.
    :math:`M = log_2(\frac{x}{y})`
    :math:`A = \frac{log_2(x * y)}{2}`

    Parameters
    ----------
    x : int or float
        Read count/density in sample 1.
    y : int or float
        Read count/density in sample 2.

    Returns
    -------
    m_value : float
        Calculated M value.
    a_value : float
        Calculated A value.
    """
    m_value = log(x, 2) - log(y, 2)
    a_value = (log(x, 2) + log(y, 2)) / 2
    return m_value, a_value


def ma_to_xy(m, a):
    """Convert (M, A) value back to read counts/densities of two samples.

    Parameters
    ----------
    m : float
        M value.
    a : float
        A vlaue.

    Returns
    -------
    x : float
        Converted read count/density in sample 1.
    y : float
        Converted read count/density in sample 2.
    """
    x = 2 ** (a + m / 2)
    y = 2 ** (a - m / 2)
    return x, y


def manorm_p(x, y):
    """Calculate MAnorm P value with given read counts/densities.

    Parameters
    ----------
    x : int or float
        Read count/density in sample 1.
    y : int or float
        Read count/density in sample 2.

    Returns
    -------
    p : float
        The probabilty of observe (`x`, `y`) given the sum `x` + `y`.
        Please refer to the manuscript of MAnorm for more information.
    """

    def _log_factorial(n):
        num = 0
        for i in range(1, n + 1):
            num += log(i)
        return num

    if x < 0 or y < 0:
        raise ValueError(f"expect x, y >= 0, got: x={x} y={y}")
    x = max(int(round(x)), 1)
    y = max(int(round(y)), 1)
    # use the log-transform to calculate p
    log_p = _log_factorial(x + y) - _log_factorial(x) - _log_factorial(y) - (
                x + y + 1) * log(2)
    if log_p < -500:
        log_p = -500
    p = exp(log_p)
    return p
