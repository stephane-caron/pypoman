#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2018 Stéphane Caron

"""Exceptions raised by pypoman.

All exceptions derive from :class:`pypoman.exceptions.Exception`, so that users
can catch any pypoman-specific error with a single ``except`` clause. Errors
that used to be raised as a built-in ``ValueError`` derive from
:class:`pypoman.exceptions.ValueError`, which subclasses both, so that existing
``except ValueError`` code keeps working.
"""

import builtins

# This module intentionally shadows the built-in ``Exception`` and
# ``ValueError`` so that pypoman errors can be caught as ``pypoman.Exception``
# and ``pypoman.ValueError`` respectively.
# pylint: disable=redefined-builtin


class Exception(builtins.Exception):
    """Base class for all exceptions raised by pypoman."""


class ValueError(Exception, builtins.ValueError):
    """A pypoman error that is also a built-in :class:`ValueError`."""


class ConversionError(ValueError):
    """A polyhedron does not match the representation it is converted to."""


class EmptyPolytopeError(ValueError):
    """A polytope or polygon is empty, i.e. its constraints are infeasible."""


class SolverError(ValueError):
    """An underlying LP or QP solver failed to return a solution."""
