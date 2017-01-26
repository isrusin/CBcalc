from distutils.core import setup, Extension
from cbclib import __version__

setup(
    name="cbclib", version=__version__,
    py_modules=["cbclib.sites", "cbclib.version"],
    ext_modules=[Extension("cbclib.counts", ["cbclib/countsmodule.c",
                                             "cbclib/countscalc.c"])]
)

