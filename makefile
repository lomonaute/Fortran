# Define program names and objects:
PGM = fem
OBJ = fedata.o link1.o plane42rect.o fea.o processor.o numeth.o main.o types.o

LIBS =  -lpgplot -lX11
LIBPATH = -Wl,-L/usr/local/lib/pgplot

# Define compiler and suffixes
COMPILE = gfortran
LINK = gfortran
FFLAGS = -fmax-errors=5 -fcheck=all -g -Wall
.SUFFIXES: .f90
%.o:%.mod

# Compile rule
.f90.o:
	$(COMPILE) $(FFLAGS) -c $<

# Link rule:
$(PGM): $(OBJ)
	$(LINK) $(FFLAGS) -o $(PGM) $(OBJ) $(LIBPATH) $(LIBS) 
	chmod go+rx $(PGM)

# Dependencies:
fea.o:         fedata.o link1.o plane42rect.o numeth.o processor.o
fedata.o:      types.o
link1.o:       types.o
main.o:        processor.o fea.o
numeth.o:      fedata.o
plane42rect.o: types.o
processor.o:   fedata.o

# Clean-up rule:
clean:
	rm -f $(OBJ) core *.inc *.vo *.mod

