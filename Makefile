CXX =  mpicxx
CXXFLAGS = -ansi -Wall -W -g  -std=c++14 -march=znver1 ${CFLAGS} -ftree-vectorize 
OPTFLAGS = -O3 -fopenmp
RPTFLAGS = -unroll-aggressive

INCLUDES = -I./headers


LIBS = -larmadillo -lstdc++ -lfftw3 -llapack -lblas -lgfortran ${LDFLAGS}

NAMEXE  = EDUS.x

OBJS  = main_MPI.o


.PHONY: exec clean



exec: $(OBJS)
	        $(CXX) $(CXXFLAGS) $(OPTFLAGS) $(RPTFLAGS) $(LFLAGS) $(OBJS) $(INCLUDES) $(LIBS) -o $(NAMEXE)
			        @echo Linked


%.o: %.cpp
	        $(CXX) $(CXXFLAGS) $(OPTFLAGS) $(RPTFLAGS) $(LFLAGS) $(INCLUDES) -c $<
			        @echo Compiled
					        
clean: 
	        $(RM) $(OBJS) $(POBJS) $(OOBJS) $(LOBJS)
	        $(RM) $(NAMEXE) $(PNAMEXE) $(ONAMEXE) $(LNAMEXE)
			$(RM) *.optrpt
