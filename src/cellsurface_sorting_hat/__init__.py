"""cellsurface_sorting_hat - sort fungal proteins into cell surface categories."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("surface_glyco")
except PackageNotFoundError:  # a checkout that was not pip-installed under this name
    __version__ = "0+unknown"
