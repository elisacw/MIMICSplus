
#-------------------------------------------------------------------------------

PROG = ./run_script_dev
LIB = /home/elisacw/anaconda3/envs/soilenv/lib
INCS = /home/elisacw/anaconda3/envs/soilenv/include

BUILD_DIR := ./bld
SRC_DIRS := ./src
MODDIR := ./mod

ifneq ($(BUILD_DIR),)
  $(shell test -d $(BUILD_DIR) || mkdir -p $(BUILD_DIR))
  FCFLAGS+= -J $(BUILD_DIR)
endif

ifneq ($(MODDIR),)
  $(shell test -d $(MODDIR) || mkdir -p $(MODDIR))
  FCFLAGS+= -J $(MODDIR)
endif

OBJ = $(BUILD_DIR)/main.o $(BUILD_DIR)/mycmimMod.o $(BUILD_DIR)/initMod.o \
	$(BUILD_DIR)/testMod.o $(BUILD_DIR)/writeMod.o $(BUILD_DIR)/paramMod.o $(BUILD_DIR)/dispmodule.o \
	$(BUILD_DIR)/shr_kind_mod.o $(BUILD_DIR)/readMod.o
#-------------------------------------------------------------------------------
FC = gfortran
FFLAGS = -g -fcheck=all -fbacktrace  -ffpe-trap=zero,invalid,overflow,underflow -O0 -Wall -ffree-line-length-'none'\
-Wl,-rpath=$(LIB) -J $(MODDIR)

#EXTRAFFLAGS =  -Wno-unused-dummy-argument -Wno-unused-parameter -Wextra -ffpe-trap=zero,invalid,overflow,underflow
FINCLUDES = -I$(INCS)

$(PROG): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) \
	-L$(LIB) -lnetcdff -lnetcdf -lmfhdf -ldf -lhdf5_hl -lhdf5 -lrt -lpthread -lz -ldl -lm -lcurl \
	-o $(PROG)

$(BUILD_DIR)/dispmodule.o: $(SRC_DIRS)/dispmodule.f90
	$(FC) $(FFLAGS) -c $(FINCLUDES) $(SRC_DIRS)/dispmodule.f90 -o $@

$(BUILD_DIR)/mycmimMod.o: $(SRC_DIRS)/mycmimMod.f90 $(BUILD_DIR)/paramMod.o $(BUILD_DIR)/dispmodule.o $(BUILD_DIR)/initMod.o $(BUILD_DIR)/writeMod.o \
$(BUILD_DIR)/testMod.o $(BUILD_DIR)/readMod.o
	$(FC) $(FFLAGS) -c $(FINCLUDES) $(SRC_DIRS)/mycmimMod.f90 -o $@

$(BUILD_DIR)/initMod.o: $(SRC_DIRS)/initMod.f90
	$(FC) $(FFLAGS) -c $(FINCLUDES) $(SRC_DIRS)/initMod.f90 -o $@

$(BUILD_DIR)/testMod.o: $(SRC_DIRS)/testMod.f90 $(BUILD_DIR)/paramMod.o $(BUILD_DIR)/dispmodule.o $(BUILD_DIR)/initMod.o
	$(FC) $(FFLAGS) -c $(FINCLUDES) $(SRC_DIRS)/testMod.f90 -o $@

$(BUILD_DIR)/writeMod.o: $(SRC_DIRS)/writeMod.f90 $(BUILD_DIR)/paramMod.o $(BUILD_DIR)/dispmodule.o $(BUILD_DIR)/initMod.o
	$(FC) $(FFLAGS) -c $(FINCLUDES) $(SRC_DIRS)/writeMod.f90 -o $@

$(BUILD_DIR)/main.o: $(SRC_DIRS)/main.f90 $(BUILD_DIR)/mycmimMod.o $(BUILD_DIR)/initMod.o
	 $(FC) $(FFLAGS) -c $(FINCLUDES) $(SRC_DIRS)/main.f90 -o $@

$(BUILD_DIR)/paramMod.o: $(SRC_DIRS)/paramMod.f90 $(BUILD_DIR)/shr_kind_mod.o $(BUILD_DIR)/initMod.o $(BUILD_DIR)/dispmodule.o
	$(FC) $(FFLAGS) -c $(FINCLUDES) $(SRC_DIRS)/paramMod.f90 -o $@

$(BUILD_DIR)/readMod.o: $(SRC_DIRS)/readMod.f90 $(BUILD_DIR)/dispmodule.o $(BUILD_DIR)/shr_kind_mod.o $(BUILD_DIR)/initMod.o
	$(FC) $(FFLAGS) -c $(FINCLUDES) $(SRC_DIRS)/readMod.f90 -o $@

$(BUILD_DIR)/shr_kind_mod.o: $(SRC_DIRS)/shr_kind_mod.f90
	$(FC) $(FFLAGS) -c $(FINCLUDES) $(SRC_DIRS)/shr_kind_mod.f90 -o $@
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
strip:


#-------------------------------------------------------------------------------
clean: strip
	rm -f $(PROG) $(BUILD_DIR)/* $(MODDIR)/*
#-------------------------------------------------------------------------------
