AUTOMAKE_OPTIONS = foreign subdir-objects
bin_PROGRAMS = medys 
medys_SOURCES =   src/all_v001/basics.f90 src/all_v001/variables.f90 src/all_v001/read_input.f90 src/all_v001/math.f90 src/all_v001/champ1.f90 src/all_v001/grille.f90 src/all_v001/NewSubroutines.f90 src/all_v001/Tests.f90 src/all_v001/string.f90 src/io_interface/io_module.f90 src/all_v001/precision.f90 src/core/error.f90 src/all_v001/matrix_modifications.f90 src/all_v001/matrix_temp_lowlvl.f90 src/all_v001/matrix_temp_highlvl.f90 src/all_v001/read_columbus.f90 src/all_v001/ortho.f90 src/all_v001/qp_integrals.f90 src/all_v001/observable.f90 src/all_v001/nouveaudrt.f90 src/all_v001/read_gamess.f90 src/all_v001/dynamique.f90  src/main/main.f90


#SUBDIRS = src/main src/all_v001
#medys_LDADD = libs/basics.o
#medys_LIBRARIES = libs/basics.o
#medys_LIBRARIES_SOURCES = src/all_v001/basics.f90
#medys_LDFLAGS = -mod include    lib/basics.o lib/variables.o lib/read_input.o lib/math.o         lib/dynamique.o  lib/read_gamess.o   lib/error.o lib/io_module.o   
#medys_LDFLAGS = -g -nogen-interface -fpe0  -traceback -debug extended -heap-arrays -fp-stack-check -debug all -openmp -warn -L$(MKLROOT)/lib/intel64 -I$(MKLROOT)/include/intel64/lp64/ -lmkl_intel_lp64 -lmkl_core  -lmkl_intel_thread -lpthread -liomp5 -lm -lmkl_blas95_lp64 -lmkl_lapack95_lp64 lib/*
#${medys_LIBRARIES}
#LIBS = -L${MKLROOT}/lib/intel64 -I${MKLROOT}/include/intel64/lp64/ -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -lpthread 
#medys_LDADD = dynamique.o
#dynamique_SOURCES = src/all_v001/dynamique.f90  src/all_v001/matrix_modifications.f90 src/all_v001/matrix_temp_highlvl.f90 src/all_v001/matrix_temp_lowlvl.f90 src/all_v001/string.f90 src/all_v001/precision.f90 src/all_v001/read_columbus.f90 src/all_v001/observable.f90 src/all_v001/nouveaudrt.f90 src/all_v001/ortho.f90 src/all_v001/qp_integrals.f90 src/all_v001/champ1.f90 src/all_v001/grille.f90 src/all_v001/NewSubroutines.f90 src/all_v001/Tests.f90 src/all_v001/basics.f90   src/all_v001/variables.f90 src/all_v001/read_input.f90 src/all_v001/math.f90
#dynamique_FCFLAGS = -mod include          lib/read_gamess.o   lib/error.o lib/io_module.o    -nogen-interface -fpe0  -traceback -debug extended -heap-arrays -fp-stack-check -debug all -openmp  -L/opt/intel/composer_xe_2013_sp1.3.174/mkl/lib/intel64 -I/opt/intel/composer_xe_2013_sp1.3.174/mkl/include/intel64/lp64/ -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -liomp5 -lm -lmkl_blas95_lp64 -lmkl_lapack95_lp64

#medys_FCFLAGS = -I${MKLROOT}/include/intel64/lp64/ -nogen-interface -fpe0 -O2 -qopenmp
#For debugging
medys_FCFLAGS = -I${MKLROOT}/include/intel64/lp64/ -nogen-interface -fpe0 -g -traceback -debug extended -heap-arrays -fp-stack-check -debug all -qopenmp
#For production
#medys_FCFLAGS = -I${MKLROOT}/include/intel64/lp64/ -O2 -nogen-interface -fpe0 -heap-arrays -fp-stack-check  -qopenmp


medys_LDADD = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -liomp5 -lmkl_blas95_lp64 -lmkl_lapack95_lp64

