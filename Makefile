#
# Makefile for building pqraaa
#

### OPTIONS ###
# These can be overridden by defining as command-line args
# Compiler type option: g++ [gnu] / icx [intel] / clang++ [amd]
CC=gnu
# Compile mode option: [debug]/[optimize]
CM=optimize
# Parallelization option: [parallel]/[serial]
#	Note: for now, serial does not actually work most of the time
CP=parallel

############################################################
#
# Compiler, compilation flags, link flags, ...
#

### COMPILER ###
# List of compilers
CCOMP_gnu = g++
CCOMP_intel = icx
CCOMP_amd = clang++
# Choose compiler with compiler option
CCOMP = $(CCOMP_$(CC))

### GLOBAL COMPILE FLAGS ###
# General compile flags
CFLAG_dep = -MMD # Tells compiler to generate dependency files (.d, located in build directory)
CFLAGS_GEN = -std=c++17 $(CFLAG_dep)
# Additional compile flags
# 
CFLAG_pqraaa_gnu   = -DPQRAAA_EIGEN_USE_BLAS -DPQRAAA_EIGEN_USE_LAPACKE -DPQRAAA_USE_ALT_LAPACK_BINDINGS
CFLAG_pqraaa_intel = -DPQRAAA_USE_MKL
CFLAG_pqraaa_amd   = -DPQRAAA_EIGEN_USE_BLAS -DPQRAAA_EIGEN_USE_LAPACKE
CFLAG_Eigen = -DEIGEN_NO_AUTOMATIC_RESIZING -DEIGEN_MPL2_ONLY
CFLAGS_ADD = $(CFLAG_pqraaa_$(CC)) $(CFLAG_Eigen)

### COMPILER-SPECIFIC FLAGS
# GNU-specific flags
CFLAGS_gnu_gen = -Wall -Wextra -march=native
CFLAGS_gnu_debug 	= -g -Og
CFLAGS_gnu_optimize = -O3 -DNDEBUG -fno-trapping-math -fno-signed-zeros -fassociative-math
CFLAGS_gnu_parallel = -fopenmp
CFLAGS_gnu_serial   = -lgomp
# Choose flags with mode and parallelization options
CFLAGS_gnu = $(CFLAGS_gnu_gen) $(CFLAGS_gnu_$(CM)) $(CFLAGS_gnu_$(CP)) 
# Intel-specific flags
CFLAGS_intel_gen = -Wall -Wremarks -w3 -xHost -fp-model precise -qmkl=sequential
CFLAGS_intel_debug 	  = -g -Og
CFLAGS_intel_optimize = -O3 -DNDEBUG
CFLAGS_intel_parallel = -qopenmp
CFLAGS_intel_serial   = -qopenmp-stubs
# Choose flags with mode and parallelization options
CFLAGS_intel = $(CFLAGS_intel_gen) $(CFLAGS_intel_$(CM)) $(CFLAGS_intel_$(CP))
# ADM-specific flags
CFLAGS_amd_gen = -Wall -Wextra -march=native
CFLAGS_amd_debug 	= -g -Og
CFLAGS_amd_optimize = -O3 -DNDEBUG
CFLAGS_amd_parallel = -fopenmp
CFLAGS_amd_serial   = -lgomp
# Choose flags with mode and parallelization options
CFLAGS_amd = $(CFLAGS_amd_gen) $(CFLAGS_amd_$(CM)) $(CFLAGS_amd_$(CP)) 
# Choose flags with compiler option
CFLAGS = $(CFLAGS_GEN) $(CFLAGS_ADD) $(CFLAGS_$(CC))

### INCLUDE FLAGS ###
# Include flags by library
#	'-isystem' is an alternative to '-I' that suppresses warnings for those included headers
IFLAG_glas2 = -isystem pqraaa-support/glas2
IFLAG_Boost = -isystem pqraaa-support/boost
IFLAG_cork  = -isystem pqraaa-support/cork
IFLAG_MUMPS =
IFLAG_Eigen = -isystem /usr/include/eigen3
# Compose all include flags
IFLAGS = -I include $(IFLAG_glas2) $(IFLAG_Boost) $(IFLAG_cork) $(IFLAG_MUMPS) $(IFLAG_Eigen)

### LINK FLAGS ###
# Link flags by library
LFLAGS_Boost  = -lboost_system -lboost_filesystem
LFLAGS_lapack_gnu   = -llapack -lblas -llapacke
LFLAGS_lapack_intel = -llapack -lblas
LFLAGS_lapack_amd   = -llapack -lblas -llapacke
# Link flags by compiler
LFLAGS_gnu = 
LFLAGS_intel = -lstdc++
LFLAGS_amd = 
# Compose all link flags
LFLAGS = $(LFLAGS_$(CC)) $(LFLAGS_Boost) $(LFLAGS_lapack_$(CC))

############################################################
#
# Define test flags
#

### TEST FLAGS ###
CFLAGS_unit 	= $(CFLAGS_GEN) $(CFLAGS_ADD) $(CFLAGS_$(CC)_gen) $(CFLAGS_$(CC)_debug) $(CFLAGS_$(CC)_parallel)
CFLAGS_perf 	= $(CFLAGS_GEN) $(CFLAGS_ADD) $(CFLAGS_$(CC)_gen) $(CFLAGS_$(CC)_optimize) $(CFLAGS_$(CC)_parallel)
CFLAGS_validate = $(CFLAGS_GEN) $(CFLAGS_ADD) $(CFLAGS_$(CC)_gen) $(CFLAGS_$(CC)_debug) $(CFLAGS_$(CC)_parallel)

############################################################
#
# Define directories
#

# C++ scripts
DIR_scripts =
# Object files 
DIR_build 	= build/
# Executables
DIR_bin 	= bin/
# Alias directory for running scripts
DIR_run 	= run/

############################################################
#
# Define pretty output (http://www.lunderberg.com/2015/08/25/cpp-makefile-pretty-output/)
#

# Colors
COLOR_BLUE   = \033[0;34m
COLOR_CYAN   = \033[0;36m
NO_COLOR    = \033[m
# Strings
COMP_STRING   = "Compile"
LINK_STRING   = "Link"
RUN_STRING    = "Run"

############################################################
#
# Define build recipes
#

clean: clean_build clean_bin


### INITIALIZATION RECIPES ###
# Create build and bin directories
init:
	@mkdir -p $(DIR_build)/examples
	@mkdir -p $(DIR_bin)/examples
	@printf "%b" "$(COLOR_BLUE)build and bin directories created.$(NO_COLOR)\n";

### INDIVIDUAL RECIPES ###
# General examples recipe
examples/%: $(DIR_build)examples/%.o
	@mkdir -p $(DIR_bin)$(dir $@)
	@printf "%b" "$(COLOR_BLUE)$(LINK_STRING)\t$(COLOR_CYAN)$(notdir $<)$(NO_COLOR)\n";
	@$(CCOMP) -o $(DIR_bin)$@ $(CFLAGS) $^ $(LFLAGS)
# General build recipe
$(DIR_build)%.o: $(DIR_scripts)%.cc
	@mkdir -p $(dir $@)
	@printf "%b" "$(COLOR_BLUE)$(COMP_STRING)\t$(COLOR_CYAN)$(notdir $<)$(NO_COLOR)\n";
	@$(CCOMP) -c -o $@ $(CFLAGS) $(IFLAGS) $<


### CLEANING RECIPES ###
clean_build:
	find $(DIR_build) -type f -delete

clean_bin:
	find $(DIR_bin) -type f -delete

### INTERMEDIATE RECIPES ###
# Line below blocks deletion of the listed intermediate files in an implicit chain rule
.SECONDARY:

### PHONY TARGETS ###
# make assumes 'FORCE' is a file that is always out of date, any recipe that depends on it will always be run
FORCE:

### INCLUDES ###
# Include dependency files
-include $(shell find $(DIR_build) -name '*.d')
