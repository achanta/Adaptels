
from cffi import FFI

def compile_adaptels_lib():
	ffibuilder = FFI()

	ffibuilder.cdef(open("adaptels.h").read())

	ffibuilder.set_source(
		"_adaptels",
		r"""
			#include "adaptels.h"
		""",
		sources=["adaptels.c"],
		library_dirs=['.'],
		extra_compile_args=['-O3', '-march=native', '-ffast-math']
		)

	ffibuilder.compile(verbose=True)


compile_adaptels_lib()