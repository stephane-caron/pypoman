# SPDX-License-Identifier: GPL-3.0-or-later

"""Other utility functions."""

from types import ModuleType

import numpy as np


def import_cdd() -> ModuleType:
    """Import pycddlib, which is an optional dependency of pypoman.

    Returns
    -------
    :
        The ``cdd`` module.

    Raises
    ------
    ImportError
        If pycddlib is not installed.
    """
    try:
        import cdd  # pylint: disable=import-outside-toplevel
    except ImportError as e:
        raise ImportError(
            "This function requires pycddlib: pip install pypoman[cdd]"
        ) from e
    return cdd


def norm(v: np.ndarray) -> float:
    """Euclidean norm.

    Parameters
    ----------
    v :
        Vector.

    Returns
    -------
    :
        Euclidean norm of the input vector.

    Notes
    -----
    This straightforward function is 2x faster than :func:`numpy.linalg.norm`
    on my machine.
    """
    return np.sqrt(np.dot(v, v))
