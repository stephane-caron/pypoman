# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added

- CICD: Documentation workflow

### Removed

- Remove cvxopt reference from polyhedron submodule

## [1.2.0] - 2025-01-06

### Changed

- Bump minimum supported Python version to 3.9
- Update pycddlib to version 3.0.0

### Fixed

- CICD: Re-enable macOS testing

### Removed

- CICD: Dropped Windows support
- Remove cdd reference from polyhedron submodule

## [1.1.1] - 2024-12-20

### Changed

- Use f-strings to shorten polygon assertion messages
- Pin pycddlib dependency version to < 3

### [1.1.0] - 2023-05-08

### Added

- Cover all submodules with unit tests

### Changed

- Raise ``ValueError`` exception when polyhedron is empty
- Remove default QP solver when projecting a point to a polytope
- Rename main branch from ``master`` to ``main``

### Fixed

- CICD: Install missing Linux dependency
- CICD: Install missing macOS dependency

### Removed

- CICD: Drop macOS runners as upstream dependencies don't support that platform

## [1.0.0] - 2023-05-18

### Added

- Add type hints to function prototypes
- Bretl: Raise an error when trying to expand a successor-less vertex
- Continuous integration for Linux, macOS and Windows
- Document all public modules
- Point to polytope projection (thanks to @peekxc)
- Unit test fixtures for vertex and halfspace enumeration

### Changed

- Convert ``setup.py`` to ``pyproject.toml``
- Drop support for Python 3.7
- Figure axes are not resized by default any more when plotting a polygon
- Remove repository-wide ``__init__.py``
- Use ``pylab.show()`` rather than IPython in examples

[unreleased]: https://github.com/qpsolvers/qpsolvers/compare/v1.2.0...HEAD
[1.2.0]: https://github.com/qpsolvers/qpsolvers/compare/v1.1.1...v1.2.0
[1.1.1]: https://github.com/qpsolvers/qpsolvers/compare/v1.1.0...v1.1.1
[1.1.0]: https://github.com/qpsolvers/qpsolvers/compare/v1.0.0...v1.1.0
[1.0.0]: https://github.com/qpsolvers/qpsolvers/releases/tag/v1.0.0
