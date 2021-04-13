# -*- coding: utf-8 -*-

"""The source code refer from "spatial_efd" module."""

__author__ = "Yuan Chang"
__copyright__ = "Copyright (C) 2016-2021"
__license__ = "AGPL"
__email__ = "pyslvs@gmail.com"

from typing import Tuple, Sized, Sequence, Union
from math import pi, sin, cos, atan2, degrees, radians
from numpy import (
    sqrt, abs, cos as np_cos, sin as np_sin, dot, sum as np_sum, array,
    ndarray, linspace, zeros, ones, diff, concatenate, cumsum,
)

_Path = Union[Sequence[Tuple[float, float]], ndarray]


def efd_fitting(path: _Path, n: int = 0) -> ndarray:
    """Curve fitting using Elliptical Fourier Descriptor.

    The path `path` will be translated to Fourier descriptor coefficients,
    then regenerate a new path as a `n` x 4 NumPy array.
    """
    contour = array(path, dtype=float)
    if n < 3:
        n = len(contour)
    harmonic = fourier_power(
        calculate_efd(contour, _nyquist(contour)),
        _nyquist(contour)
    )
    coeffs = calculate_efd(contour, harmonic)
    coeffs, rotation = normalize_efd(coeffs, size_invariant=False)
    locus_v = locus(contour)
    # New path
    contour = inverse_transform(coeffs, locus_v, n, harmonic)
    return rotate_contour(contour, -rotation, locus_v)


def normalize_efd(
    coeffs: ndarray,
    size_invariant: bool = True
) -> Tuple[ndarray, float]:
    """
    Normalize the Elliptical Fourier Descriptor coefficients for a polygon.

    Implements Kuhl and Giardina method of normalizing the coefficients
    An, Bn, Cn, Dn. Performs 3 separate normalizations. First, it makes the
    data location invariant by re-scaling the data to a common origin.
    Secondly, the data is rotated with respect to the major axis. Thirdly,
    the coefficients are normalized with regard to the absolute value of A_1.
    This code is adapted from the pyefd module. See the original paper for
    more detail:

    Kuhl, FP and Giardina, CR (1982). Elliptic Fourier features of a closed
    contour. Computer graphics and image procesnp_sing, 18(3), 236-258.

    Args:
        coeffs: A numpy array of shape (n, 4) representing the
            four coefficients for each harmonic computed.
        size_invariant: Set to True (the default) to perform the third
            normalization and false to return the data withot this procesnp_sing
            step. Set this to False when plotting a comparison between the
            input data and the Fourier ellipse.
    Returns:
        A tuple consisting of a numpy.ndarray of shape (harmonics, 4)
        representing the four coefficients for each harmonic computed and
        the rotation in degrees applied to the normalized contour.
    """
    # Make the coefficients have a zero phase shift from
    # the first major axis. Theta_1 is that shift angle.
    theta1 = 0.5 * atan2(
        2 * (coeffs[0, 0] * coeffs[0, 1] + coeffs[0, 2] * coeffs[0, 3]),
        coeffs[0, 0] ** 2 - coeffs[0, 1] ** 2
        + coeffs[0, 2] ** 2 - coeffs[0, 3] ** 2
    )
    # Rotate all coefficients by theta_1
    for n in range(coeffs.shape[0]):
        angle = (n + 1) * theta1
        coeffs[n, :] = dot(
            array([
                [coeffs[n, 0], coeffs[n, 1]],
                [coeffs[n, 2], coeffs[n, 3]],
            ]),
            array([
                [np_cos(angle), -np_sin(angle)],
                [np_sin(angle), np_cos(angle)],
            ])
        ).flatten()
    # Make the coefficients rotation invariant by rotating so that
    # the semi-major axis is parallel to the x-axis.
    psi1 = atan2(coeffs[0, 2], coeffs[0, 0])
    psi2 = array([
        [np_cos(psi1), np_sin(psi1)],
        [-np_sin(psi1), np_cos(psi1)],
    ])
    # Rotate all coefficients by -psi_1.
    for n in range(coeffs.shape[0]):
        rot = array([
            [coeffs[n, 0], coeffs[n, 1]],
            [coeffs[n, 2], coeffs[n, 3]],
        ])
        coeffs[n, :] = psi2.dot(rot).flatten()
    if size_invariant:
        # Obtain size-invariance by normalizing.
        coeffs /= abs(coeffs[0, 0])
    return coeffs, degrees(psi1)


def locus(contour: ndarray) -> Tuple[float, float]:
    """
    Compute the dc coefficients, used as the locus when calling
    inverse_transform().

    This code is adapted from the pyefd module. See the original paper for
    more detail:

    Kuhl, FP and Giardina, CR (1982). Elliptic Fourier features of a closed
    contour. Computer graphics and image procesnp_sing, 18(3), 236-258.

    Args:
        contour: A n x 2 numpy array represents a path.
    Returns:
        A tuple containing the c and d coefficients.
    """
    dxy = diff(contour, axis=0)
    dt = sqrt((dxy ** 2).sum(axis=1))
    t = concatenate(([0], cumsum(dt)))
    zt = t[-1]
    diffs = diff(t ** 2)
    xi = cumsum(dxy[:, 0]) - dxy[:, 0] / dt * t[1:]
    a0 = np_sum(dxy[:, 0] / (2 * dt) * diffs + xi * dt) / (zt + 1e-20)
    delta = cumsum(dxy[:, 1]) - dxy[:, 1] / dt * t[1:]
    c0 = np_sum(dxy[:, 1] / (2 * dt) * diffs + delta * dt) / (zt + 1e-20)
    # A0 and CO relate to the first point of the contour array as origin
    # Adding those values to the coefficients to make them relate to true origin
    return contour[0, 0] + a0, contour[0, 1] + c0


def calculate_efd(contour: ndarray, harmonic: int = 10) -> ndarray:
    """
    Compute the Elliptical Fourier Descriptors for a polygon.

    Implements Kuhl and Giardina method of computing the coefficients
    An, Bn, Cn, Dn for a specified number of harmonics. This code is adapted
    from the pyefd module. See the original paper for more detail:

    Kuhl, FP and Giardina, CR (1982). Elliptic Fourier features of a closed
    contour. Computer graphics and image procesnp_sing, 18(3), 236-258.

    Args:
        contour: A n x 2 numpy array represents a path.
        harmonic: The number of harmonics to compute for the given
            shape, defaults to 10.
    Returns:
        A numpy array of shape (harmonics, 4) representing the
        four coefficients for each harmonic computed.
    """
    dxy = diff(contour, axis=0)
    dt = sqrt((dxy ** 2).sum(axis=1))
    t = concatenate(([0], cumsum(dt)))
    zt = t[-1]
    phi = (2. * pi * t) / (zt + 1e-20)
    coeffs = zeros((harmonic, 4))
    for n in range(1, harmonic + 1):
        const = zt / (2 * n * n * pi * pi)
        phi_n = phi * n
        d_np_cos_phi_n = np_cos(phi_n[1:]) - np_cos(phi_n[:-1])
        d_np_sin_phi_n = np_sin(phi_n[1:]) - np_sin(phi_n[:-1])
        coeffs[n - 1, :] = (
            const * np_sum(dxy[:, 1] / dt * d_np_cos_phi_n),
            const * np_sum(dxy[:, 1] / dt * d_np_sin_phi_n),
            const * np_sum(dxy[:, 0] / dt * d_np_cos_phi_n),
            const * np_sum(dxy[:, 0] / dt * d_np_sin_phi_n),
        )
    return coeffs


def inverse_transform(
    coeffs: ndarray,
    locus_v: Tuple[float, float] = (0., 0.),
    n: int = 300,
    harmonic: int = 10
) -> ndarray:
    """
    Perform an inverse fourier transform to convert the coefficients back into
    spatial coordinates.

    Implements Kuhl and Giardina method of computing the performing the
    transform for a specified number of harmonics. This code is adapted
    from the pyefd module. See the original paper for more detail:

    Kuhl, FP and Giardina, CR (1982). Elliptic Fourier features of a closed
    contour. Computer graphics and image procesnp_sing, 18(3), 236-258.

    Args:
        coeffs: A numpy array of shape (n, 4) representing the
            four coefficients for each harmonic computed.
        locus_v: The x,y coordinates of the centroid of the contour being
            generated. Use calculate_dc_coefficients() to generate the correct
            locus for a shape.
        n: The number of coordinate pairs to compute. A larger value will
            result in a more complex shape at the expense of increased
            computational time. Defaults to 300.
        harmonic: The number of harmonics to be used to generate
            coordinates, defaults to 10. Must be <= coeffs.shape[0]. Supply a
            smaller value to produce coordinates for a more generalized shape.
    Returns:
        A n x 2 numpy array represents a contour.
    """
    t = linspace(0, 1, n)
    contour = ones((n, 2), dtype=float)
    contour[:, 0] *= locus_v[0]
    contour[:, 1] *= locus_v[1]
    for n in range(harmonic):
        angle = 2 * (n + 1) * pi * t
        cosine = np_cos(angle)
        sine = np_sin(angle)
        contour[:, 0] += coeffs[n, 2] * cosine + coeffs[n, 3] * sine
        contour[:, 1] += coeffs[n, 0] * cosine + coeffs[n, 1] * sine
    return contour


def _nyquist(zx: Sized) -> int:
    """
    Returns the maximum number of harmonics that can be computed for a given
    contour, the Nyquist Frequency.

    See this paper for details:
    C. np_costa et al. / Postharvest Biology and Technology 54 (2009) 38-47

    Args:
        zx (list): A list (or numpy array) of x coordinate values.
    Returns:
        int: The nyquist frequency, expressed as a number of harmonics.
    """
    return len(zx) // 2


def fourier_power(
    coeffs: ndarray,
    nyq: int,
    threshold: float = 0.9999
) -> int:
    """
    Compute the total Fourier power and find the minimum number of harmonics
    required to exceed the threshold fraction of the total power.

    This is a good method for identifying the number of harmonics to use to
    describe a polygon. For more details see:

    C. np_costa et al. / Postharvest Biology and Technology 54 (2009) 38-47

    Warning:
        The number of coefficients must be >= the Nyquist Frequency.
    Args:
        coeffs: A numpy array of shape (n, 4) representing the
            four coefficients for each harmonic computed.
        nyq: The Nyquist Frequency.
        threshold: The threshold fraction of the total Fourier power,
            the default is 0.9999.
    Returns:
        The number of harmonics required to represent the contour above
        the threshold Fourier power.
    """
    total_power = 0
    current_power = 0
    for i in range(nyq):
        total_power += 0.5 * np_sum(coeffs[i, :] ** 2)
    for i in range(nyq):
        current_power += 0.5 * np_sum(coeffs[i, :] ** 2)
        if current_power / total_power > threshold:
            return i + 1
    return nyq


def rotate_contour(
    contour: ndarray,
    rotation: float,
    centroid: Tuple[float, float]
) -> ndarray:
    """
    Rotates a contour about a point by a given amount expressed in degrees.

    Operates by calling rotatePoint() on each x,y pair in turn. X and Y must
    have the same dimensions.

    Args:
        contour: A n x 2 numpy array represents a path.
        rotation: The angle in degrees for the contour to be rotated by.
        centroid: A tuple containing the x,y coordinates of the centroid to
            rotate the contour about.
    Returns:
        A n x 2 numpy array represents a contour.
    """
    angle = radians(rotation)
    cpx, cpy = centroid
    out = zeros(contour.shape, dtype=contour.dtype)
    for i in range(len(contour)):
        dx = contour[i, 0] - cpx
        dy = contour[i, 1] - cpy
        out[i] = (
            cpx + dx * cos(angle) - dy * sin(angle),
            cpy + dx * sin(angle) + dy * cos(angle)
        )
    return out
