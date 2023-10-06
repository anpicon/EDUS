CXX = mpiicpc
CXXFLAGS = -ansi -Wall -W -g  -std=c++14  
OPTFLAGS = -O3 -no-prec-div  -qopenmp  
RPTFLAGS = -unroll-aggressive -ansi-alias 
LFLAGS = -L/usr/local/Intel/impi/2018.2.199/intel64/lib   -L/usr/local/fftw/lib -L/usr/local/Intel/lib/intel64/ -L/usr/local/Intel/lib/intel64/ -L/opt/intel/oneapi/mkl/2022.2.1/lib/intel64
INCLUDES = -I/lustre/local/intel/composer-2019/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin/ -I/lustre/local/intel/composer-2019/compilers_and_libraries_2018.2.199/linux/mpi/intel64/include -I./headers -I/lustre/local/fftw/include -I/usr/local/fftw/include/ -I/usr/local/armadillo/9.800.2/include
LIBS = -larmadillo -lstdc++ -lfftw3  -lmkl_rt
NAMEXE  = W.x
ONAMEXE = mpicbwe.x
PNAMEXE = Core.x
LNAMEXE = dipolecc.x
OBJS  = cBWE.o
OOBJS = main_MPI.o
POBJS = Core_2.o
LOBJS = dipole_cc.o

.PHONY: exec clean

oexec: $(OOBJS)
	        $(CXX) $(CXXFLAGS) $(OPTFLAGS) $(RPTFLAGS) $(LFLAGS) $(OOBJS) $(INCLUDES) $(LIBS) -o $(ONAMEXE)
			        @echo Linked

exec: $(OBJS)
	        $(CXX) $(CXXFLAGS) $(OPTFLAGS) $(RPTFLAGS) $(LFLAGS) $(OBJS) $(INCLUDES) $(LIBS) -o $(NAMEXE)
			        @echo Linked

pexec: $(POBJS)
	        $(CXX) $(CXXFLAGS) $(OPTFLAGS) $(RPTFLAGS) $(LFLAGS) $(POBJS) $(INCLUDES) $(LIBS) -o $(PNAMEXE)
			        @echo Linked


lexec: $(LOBJS)
	        $(CXX) $(CXXFLAGS) $(OPTFLAGS) $(RPTFLAGS) $(LFLAGS) $(LOBJS) $(INCLUDES) $(LIBS) -o $(LNAMEXE)
			        @echo Linked
%.o: %.cpp
	        $(CXX) $(CXXFLAGS) $(OPTFLAGS) $(RPTFLAGS) $(LFLAGS) $(INCLUDES) -c $<
			        @echo Compiled
					        
clean: 
	        $(RM) $(OBJS) $(POBJS) $(OOBJS) $(LOBJS)
	        $(RM) $(NAMEXE) $(PNAMEXE) $(ONAMEXE) $(LNAMEXE)
			$(RM) *.optrpt
