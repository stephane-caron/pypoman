#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2018-2020 Stephane Caron <stephane.caron@normalesup.org>
#
# This file is part of pypoman <https://github.com/stephane-caron/pypoman>.
#
# pypoman is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# pypoman is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# pypoman. If not, see <http://www.gnu.org/licenses/>.

"""Other utility functions."""

from datetime import datetime

import numpy as np


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


def normalize(v: np.ndarray) -> np.ndarray:
    """Normalize a vector.

    Parameters
    ----------
    v :
        Any vector.

    Returns
    -------
    :
        Unit vector directing `v`.

    Notes
    -----
    This method doesn't catch ``ZeroDivisionError`` exceptions on purpose.
    """
    return v / norm(v)


def warn(msg: str) -> None:
    """
    Log a warning message (in yellow) to stdout.

    Parameters
    ----------
    msg :
        Warning message.
    """
    date = datetime.now().strftime("%Y-%m-%d %H:%M:%S,%f")[:-3]
    print("%c[0;%d;48m%s pypoman [WARN] %s%c[m" % (0x1B, 33, date, msg, 0x1B))
