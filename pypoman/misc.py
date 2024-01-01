#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2018 Stéphane Caron

"""Other utility functions."""

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
