from distutils.core import setup, Extension

setup(
    name="cbclib", version="0.9",
    py_modules=["cbclib.sites"],
    ext_modules=[Extension("cbclib.counts", ["cbclib/countsmodule.c", "cbclib/countscalc.c"])]
)

