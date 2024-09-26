from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import subprocess
import os

os.environ["CC"] = "gcc"

VERSION = "1.0.0"

# Define the sources and the corresponding C flags
sources = [
    "src/booster/src/hashtables_bfields.c",
    "src/booster/src/M1M2_hashmap.c",
    "src/booster/src/tree.c",
    "src/booster/src/stats.c",
    "src/booster/src/prng.c",
    "src/booster/src/hashmap.c",
    "src/booster/src/version.c",
    "src/booster/src/sort.c",
    "src/booster/src/io.c",
    "src/booster/src/tree_utils.c",
    "src/booster/src/bitset_index.c",
    "src/booster/src/rapid_transfer.c",
    "src/booster/src/heavy_paths.c",
    "src/booster/src/debug.c",
    "src/booster/src/kludge.c",
    "src/booster/src/node_stack.c",
    "src/booster/src/split.c",
    "src/booster/src/booster.c",  # Your main file
]

# Define the compilation flags
# compile_flags = ["-Wall", "-O0", "-g", "-fopenmp", "-fsanitize=leak", "-fPIC", "-DDEBUG", f'-DVERSION="{VERSION}"']
compile_flags = ["-Wall", "-O0", "-g", "-fopenmp", "-fsanitize=leak", "-fPIC", "-DNDEBUG", f'-DVERSION="{VERSION}"']
# compile_flags = ["-Wall", "-O3", "-fopenmp", "-fPIC", "-DNDEBUG", f'-DVERSION="{VERSION}"']
link_flags = ["-lm", "-fopenmp"]

booster_extension = Extension(
    "booster",  # Name of the extension
    sources=sources,
    extra_compile_args=compile_flags,
    extra_link_args=link_flags,
)

setup(
    name="Consensus",
    version="0.0.2",
    ext_modules=[booster_extension],  # Add the extension here
)