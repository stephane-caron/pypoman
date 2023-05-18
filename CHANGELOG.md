# Changelog

All notable changes to this project will be documented in this file.

## Unreleased

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
