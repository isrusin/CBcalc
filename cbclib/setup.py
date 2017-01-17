from distutils.core import setup, Extension

setup(
    name="counts", version="0.1",
    ext_modules=[Extension("counts", ["countsmodule.c", "countscalc.c"])]
)

