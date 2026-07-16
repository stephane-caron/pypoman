#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2018 Stéphane Caron

"""Exceptions raised by pypoman.

All exceptions derive from :class:`pypoman.exceptions.Exception`, so that users
can catch any pypoman-specific error with a single ``except`` clause.
"""

import builtins

# This module intentionally shadows the built-in ``Exception`` so that pypoman
# errors can be caught as ``pypoman.Exception``.
# pylint: disable=redefined-builtin


class Exception(builtins.Exception):
    """Base class for all exceptions raised by pypoman."""


class ConversionError(Exception):
    """A polyhedron does not match the representation it is converted to."""


class EmptyPolytopeError(Exception):
    """A polytope or polygon is empty, i.e. its constraints are infeasible."""


class SolverError(Exception):
    """An underlying LP or QP solver failed to return a solution."""
