CCF = -DHDF
FC = -c
FFLAGS = -Wall -Wextra -Wconversion -Wno-maybe-uninitialized -pedantic -ffree-line-length-none -fbounds-check
FOPT = -O3 -march=native -funroll-loops -flto
FOMP = -fopenmp 

OBJECTS_DIR = src/objects
BIN_DIR = bin
EXECUTABLE = $(BIN_DIR)/run

F90_FILES := $(wildcard src/*.f90)
FCOMP = h5fc

OBJECTS := $(patsubst   src/%.f90, $(OBJECTS_DIR)/%.o, $(F90_FILES))

$(EXECUTABLE): $(OBJECTS)
		@mkdir -p $(BIN_DIR)
		$(FCOMP) -cpp  $(CCF) $(FFLAGS) $(FOMP) $(FOPT) -L/usr/local/lib/ $^ -o $@ -llapack -lblas 
		@mv *.mod bin

$(OBJECTS_DIR)/%.o : src/%.f90
		@mkdir -p $(OBJECTS_DIR)
		$(FCOMP) -cpp $(CCF) $(FC) $(FFLAGS) $(FOPT) $(FOMP) -L/usr/local/lib/ $< -o $@ -llapack -lblas 
		

clean:
		rm -rf $(OBJECTS_DIR) 
		rm -rf $(EXECUTABLE)
		rm -rf $(BIN_DIR)