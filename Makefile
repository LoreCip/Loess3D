USEMPI ?= -DUSE_MPI

FC = -c
FFLAGS = -Wall -Wextra -Wconversion -Wno-maybe-uninitialized -pedantic -ffree-line-length-none -fbounds-check
FOPT = -O3 -march=native -funroll-loops -flto
FOMP = -fopenmp 
# Check if USE_MPI is defined, and include MPI flags accordingly
ifeq ($(USEMPI), -DUSE_MPI)
    FMPI = -ffile-prefix-map=/build/mpich-0xgrG5/mpich-4.0=. -flto=auto -ffat-lto-objects -fstack-protector-strong -fallow-invalid-boz -fallow-argument-mismatch -Wl,-Bsymbolic-functions -flto=auto -ffat-lto-objects -flto=auto -Wl,-z,relro -I/usr/include/x86_64-linux-gnu/mpich -L/usr/lib/x86_64-linux-gnu -lmpichfort -lmpich
else
    FMPI =
endif

OBJECTS_DIR = src/objects
BIN_DIR = bin
EXECUTABLE = $(BIN_DIR)/run

F90_FILES := src/modules/ioH5.f90 src/modules/sort_interface.f90 src/modules/utils.f90 src/modules/math.f90 src/main.f90
FCOMP ?= h5fc

OBJECTS := $(patsubst src/modules/%.f90, $(OBJECTS_DIR)/%.o, $(F90_FILES))

$(EXECUTABLE): $(OBJECTS)
		@mkdir -p $(BIN_DIR)
		$(FCOMP) -cpp $(USEMPI) $(FFLAGS) $(FOMP) $(FMPI) $(FOPT) $^ -o $@ -llapack -lblas 
		@mv *.mod src/modules/
		@mv main.o src/objects
		
$(OBJECTS_DIR)/%.o : src/modules/%.f90
		@mkdir -p $(OBJECTS_DIR)
		$(FCOMP) -cpp $(USEMPI) $(FC) $(FFLAGS) $(FOPT) $(FMPI) $(FOMP) $< -o $@ -llapack -lblas 
		

clean:
		rm -rf $(OBJECTS_DIR) 
		rm -rf $(EXECUTABLE)
		rm -rf $(BIN_DIR)
		rm -f src/modules/*.mod