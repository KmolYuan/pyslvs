cimport cython
from numpy import zeros, array, arange, interp, float64 as np_float
from libc.math cimport cos, sin, fabs, sqrt, atan2, INFINITY as INF
from .expression cimport Coord
from .metaheuristics.utility cimport ObjFunc


def norm_path(path, scale=1):
    """Python wrapper of normalization function."""
    cdef Coord[:] path_m = array([
        Coord.__new__(Coord, x, y) for x, y in path], dtype=object)
    _norm(path_m, scale)
    return [(c.x, c.y) for c in path_m]


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _norm(Coord[:] path, double scale):
    """Normalization implementation inplace."""
    cdef Coord centre = Coord.__new__(Coord, 0, 0)
    cdef Coord c
    for c in path:
        centre.x += c.x
        centre.y += c.y
    cdef size_t end = len(path)
    centre.x /= end
    centre.y /= end
    cdef double[:] angle = zeros(end + 1, dtype=np_float)
    cdef double[:] length = zeros(end + 1, dtype=np_float)
    cdef size_t sp = 0
    cdef size_t i
    for i in range(len(path)):
        c = path[i]
        angle[i] = atan2(c.y - centre.y, c.x - centre.x)
        length[i] = centre.distance(c)
        if length[i] > length[end]:
            length[end] = length[i]
            angle[end] = angle[i]
            sp = end
    _aligned(path, sp)
    cdef double[:] bound = array([INF, -INF], dtype=np_float)
    cdef double a
    for i in range(len(path)):
        c = path[i]
        a = angle[i] - angle[end]
        c.x = length[i] * cos(a)
        c.y = length[i] * sin(a)
        if c.x < bound[0]:
            bound[0] = c.x
        if c.x > bound[1]:
            bound[1] = c.x
    scale /= (bound[1] - bound[0])
    for c in path:
        c.x *= scale
        c.y *= scale


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void _aligned(Coord[:] path, size_t sp):
    """Split 1D path from sp, concatenate to end."""
    if sp == 0:
        return
    cdef double[:, :] tmp = zeros((sp, 2), dtype=np_float)
    cdef size_t i
    for i in range(sp):
        tmp[i, 0] = path[i].x
        tmp[i, 1] = path[i].y
    for i in range(sp, len(path)):
        path[i - sp].x = path[i].x
        path[i - sp].y = path[i].y
    for i in range(sp):
        path[len(path) - sp + i].x = tmp[i, 0]
        path[len(path) - sp + i].y = tmp[i, 1]


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double _mean(double[:] p):
    """Calculate mean of the memory view."""
    cdef double s = 0
    cdef double v
    for v in p:
        s += v
    return s / len(p)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:, :] _to_numpy(Coord[:] path):
    """To memory view."""
    cdef double[:, :] p = zeros((len(path), 2), dtype=np_float)
    cdef int i
    cdef Coord c
    for i, c in enumerate(path):
        p[i, 0] = c.x
        p[i, 1] = c.y
    return p


@cython.boundscheck(False)
@cython.wraparound(False)
cdef Coord[:] _to_coord(double[:, :] path):
    """To coordinates."""
    coords = []
    cdef int i
    for i in range(len(path)):
        coords.append(Coord.__new__(Coord, path[i, 0], path[i, 1]))
    cdef Coord[:] path_m = array(coords, dtype=object)
    return path_m


def curvature(path):
    r"""Calculate the signed curvature and return as an array.

    $$
    \kappa(t) = \frac{x'y'' - x''y'}{(x'^2 + y'^2)^\frac{3}{2}}
    $$
    """
    cdef Coord[:] path_m = array([
        Coord.__new__(Coord, x, y) for x, y in path], dtype=object)
    return array(_curvature(path_m))


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:] _curvature(Coord[:] path):
    """Calculate the signed curvature."""
    cdef double[:, :] p = _to_numpy(path)
    cdef double[:, :] p1d = _derivative(p)
    cdef double[:, :] p2d = _derivative(p1d)
    cdef double[:] k = zeros(len(path) - 2, dtype=np_float)
    cdef int i
    for i in range(len(path) - 2):
        k[i] = ((p1d[i, 0] * p2d[i, 1] - p2d[i, 0] * p1d[i, 1])
                / (p1d[i, 0] * p1d[i, 0] + p1d[i, 1] * p1d[i, 1]) ** 1.5)
    return k


def derivative(double[:, :] p):
    """Differential function. Return $p'$."""
    return array(_derivative(p))


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:, :] _derivative(double[:, :] p):
    """Differential function backend."""
    cdef double[:, :] pd = zeros((len(p), 2), dtype=np_float)
    cdef double max0 = 0
    cdef double max1 = 0
    cdef int i, j
    for i in range(len(p)):
        j = i + 1
        if j >= len(p):
            j = 0
        pd[i, 0] = p[j, 0] - p[i, 0]
        pd[i, 1] = p[j, 1] - p[i, 1]
        if pd[i, 0] > max0:
            max0 = pd[i, 0]
        if pd[i, 1] > max1:
            max1 = pd[i, 1]
    i = len(p) - 1
    if max0 == pd[i, 0] or max1 == pd[i, 1]:
        return pd[:i]
    else:
        return pd


def path_signature(double[:] k):
    r"""Require a curvature, return path signature.
    It's composed by curvature $\kappa$ and a $K$ value.

    $$
    K = \int^t_0 |\kappa(t)| dt
    $$

    >>> path_signature(curvature(...))
    """
    return array(_path_signature(k))


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:, :] _path_signature(double[:] k):
    """Require a curvature, return path signature."""
    cdef double[:, :] s = zeros((len(k), 2), dtype=np_float)
    cdef int i
    for i in range(len(k)):
        if i > 0:
            s[i, 0] = s[i - 1, 0]
        s[i, 0] += fabs(k[i])
    s[:, 1] = k
    return s


def cross_correlation(double[:, :] p1, double[:, :] p2, double t):
    r"""Compare two path signature and return as an 1d array.

    $$
    \begin{aligned}
    C_n(j, W, P) &= \left|\sum_i^{l_P} \frac{(W_{i + j}
    - \overline{W}_{j\rightarrow j + l_P})(P_i-\overline{P})}{
    \sqrt{\sum_i^{l_P}(W_{i + j} - \overline{W}_{j\rightarrow j + l_P})^2
    \sum_i^{l_P}(P_i - \overline{P})^2}}\right|
    \\
    S &= \arg\max\{C_n(j)\} t
    \end{aligned}
    $$

    >>> ps1 = path_signature(curvature(...))
    >>> ps2 = path_signature(curvature(...))
    >>> cc = cross_correlation(ps1, ps2, len(ps1))
    """
    return array(_cross_correlation(p1, p2, t), dtype=np_float)


@cython.boundscheck(False)
@cython.wraparound(False)
cdef double[:] _cross_correlation(double[:, :] ps1, double[:, :] ps2, double t):
    """Compare two path signature."""
    cdef double[:] p1 = interp(arange(0, ps1[len(ps1) - 1, 0], t), ps1[:, 0],
                               ps1[:, 1])
    cdef double[:] p2 = interp(arange(0, ps2[len(ps2) - 1, 0], t), ps2[:, 0],
                               ps2[:, 1])
    cdef int diff = len(p1) - len(p2)
    cdef double[:] cn = zeros(diff, dtype=np_float)
    cdef int i, j, k
    cdef double m1, m2, tmp, tmp1, tmp2
    for j in range(diff):
        for i in range(len(p2)):
            m1 = _mean(p1[j:j + len(p2)])
            m2 = _mean(p2)
            tmp1 = 0
            for k in range(len(p2)):
                tmp = p1[k + j] - m1
                tmp1 += tmp * tmp
            tmp2 = 0
            for k in range(len(p2)):
                tmp = p2[k] - m2
                tmp2 += tmp * tmp
            cn[j] += (p1[i + j] - m1) * (p2[i] - m2) / sqrt(tmp1 * tmp2)
        cn[j] = fabs(cn[j])
    return cn

'''
@cython.final
@cython.boundscheck(False)
@cython.wraparound(False)
cdef class Curvature(ObjFunc):
    """The objective function to compare curvature."""

    def __cinit__(self, dict mech):
        # mech = {
        #     'expression': List[VPoint],
        #     'input': OrderedDict([((b0, d0), [start, end]), ...]),
        #     'placement': {pt: (x, y, r)},
        #     'target': {pt: [(x0, y0), (x1, y1), ...]},
        #     'same': {pt: match_to_pt},
        #     # Bounds of base link length
        #     'upper': float,
        #     'lower': float,
        #     'shape_only': bool,
        # }
        placement = mech.get('placement', {})
        if len(placement) == 0:
            raise ValueError("no grounded joint")
        target = mech.get('target', {})
        if len(target) == 0:
            raise ValueError("no target joint")
        check_set = {len(t) for t in target.values()}
        if len(check_set) != 1:
            raise ValueError("target paths should be in the same size")
        self.target_count = check_set.pop()
'''
