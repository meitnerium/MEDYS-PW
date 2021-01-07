testGS: testGS.f90 ortho.o
	ifort -g -nogen-interface -check all -fpe0  -traceback -debug extended -heap-arrays -fp-stack-check -debug all -openmp  -mod ../../include    ../../lib/basics.o ../../lib/variables.o ../../lib/read_input.o ../../lib/math.o ../../lib/Tests.o ../../lib/NewSubroutines.o ../../lib/grille.o ../../lib/champ1.o ../../lib/qp_integrals.o ../../lib/ortho.o ../../lib/nouveaudrt.o ../../lib/observable.o ../../lib/dynamique.o ../../lib/read_columbus.o     ../../lib/precision.o ../../lib/string.o ../../lib/error.o     ../../lib/io_module.o     ../../lib/matrix_temp_lowlvl.o ../../lib/matrix_temp_highlvl.o ../../lib/matrix_modifications.o       testGS.f90 -o testGS -L/opt/intel/composer_xe_2013_sp1.3.174/mkl/../../lib/intel64 -I/opt/intel/composer_xe_2013_sp1.3.174/mkl/include/intel64/lp64/ -lmkl_intel_lp64 -lmkl_core  -lmkl_intel_thread -lpthread -liomp5 -lm -lmkl_blas95_lp64 -lmkl_lapack95_lp64
all:	testGS
clean:
	rm testGS ortho.mod
	#ifort -mod ../../include ../../lib/ortho.o testGS.f90 -o testGS
