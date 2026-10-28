# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added

- Support for Python 3.15

### Changed

- CICD: Run macOS tests on `macos-latest`
- CICD: Update actions-gh-pages to v4
- CICD: Update checkout actions to v7
- CICD: Update setup-micromamba actions to v3,
- Tooling: install cvxopt from conda-forge in pixi environments
- Tooling: move pixi manifest from `pixi.toml` to `pyproject.toml`

### Fixed

- Tooling: remove leftover `pixi.toml`

### Removed

- Drop support for Python 3.10 (reached end of life)

## [2.0.0] - 2026-07-17

### Added

- CICD: Documentation workflow
- Example: Chebyshev center of a polygon
- Exceptions: ``ConversionError``, ``EmptyPolytopeError`` and ``SolverError``
- Exceptions: ``pypoman.Exception`` base class
- Support for Python 3.13 and 3.14
- Tooling: pixi manifest to run the CI locally (`pixi run tests`, `pixi run lint`)

### Changed

- CICD: Migrate development tooling from tox to pixi
- CICD: Update checkout action to v4
- CICD: Use ruff as code formatter in place of black
- **Breaking:** All library errors now derive from ``pypoman.Exception``

### Removed

- Drop support for Python 3.8 and 3.9 (reached end of life)
- Remove cvxopt reference from polyhedron submodule

### Fixed

- **Breaking:** Raise an error when point-to-polytope projection finds no solution, rather than returning ``None``

## [1.2.0] - 2025-01-06

### Changed

- Bump minimum supported Python version to 3.9
- Update pycddlib to version 3.0.0

### Removed

- CICD: Dropped Windows support
- Remove cdd reference from polyhedron submodule

### Fixed

- CICD: Re-enable macOS testing

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

### Removed

- CICD: Drop macOS runners as upstream dependencies don't support that platform

### Fixed

- CICD: Install missing Linux dependency
- CICD: Install missing macOS dependency

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

[unreleased]: https://github.com/qpsolvers/qpsolvers/compare/v2.0.0...HEAD
[2.0.0]: https://github.com/qpsolvers/qpsolvers/releases/tag/v2.0.0
[1.2.0]: https://github.com/qpsolvers/qpsolvers/releases/tag/v1.2.0
[1.1.1]: https://github.com/qpsolvers/qpsolvers/releases/tag/v1.1.1
[1.1.0]: https://github.com/qpsolvers/qpsolvers/releases/tag/v1.1.0
[1.0.0]: https://github.com/qpsolvers/qpsolvers/releases/tag/v1.0.0
