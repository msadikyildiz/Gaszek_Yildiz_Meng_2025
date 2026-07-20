from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("fitness-landscape-graph")
except PackageNotFoundError:
    # Package not installed (rare in your workflow, but good to have a fallback)
    __version__ = "0.0.0"

__all__ = [
    "__version__",
    # add other public symbols here later, e.g. "SomeModel", "train", ...
]
