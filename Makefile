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
#CFLAG_BEACHpack_gnu   = -DBEACHPACK_EIGEN_USE_BLAS -DBEACHPACK_EIGEN_USE_LAPACKE -DBEACHPACK_USE_ALT_LAPACK_BINDINGS
#CFLAG_BEACHpack_intel = -DBEACHPACK_USE_MKL
#CFLAG_BEACHpack_amd   = -DBEACHPACK_EIGEN_USE_BLAS -DBEACHPACK_EIGEN_USE_LAPACKE
#CFLAG_Eigen = -DEIGEN_NO_AUTOMATIC_RESIZING -DEIGEN_MPL2_ONLY
#CFLAGS_ADD = $(CFLAG_BEACHpack_$(CC)) $(CFLAG_Eigen)

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

### BUILD RECIPE SETS ###
all: experiments/unfcSphereAnalysis

# Clean build and binary files
clean: clean_build clean_bin

# Gather all test source files and executables
SRC_unit 	 = $(wildcard test/unit/*.cc)
EXE_unit 	 = $(patsubst %.cc,$(DIR_run)%,$(SRC_unit))
SRC_perf 	 = $(wildcard test/perf/*.cc)
EXE_perf 	 = $(patsubst %.cc,$(DIR_run)%,$(SRC_perf))
SRC_validate = $(wildcard test/validate/*.cc)
EXE_validate = $(patsubst %.cc,$(DIR_run)%,$(SRC_validate))

# Build and run unit tests
unit: $(EXE_unit)
# unit: 	$(shell find $(DIR_scripts)test/unit -name '*.cc' -exec basename {} \; | sed 's/.cc//g' | sed 's/^/test\/unit\//g')

# Build and run performance tests
perf: $(EXE_perf)
# perf: 	$(shell find $(DIR_scripts)test/perf -name '*.cc' -exec basename {} \; | sed 's/.cc//g' | sed 's/^/test\/perf\//g')

# Build and run validation tests
validate: $(EXE_validate)
# validate: $(shell find $(DIR_scripts)test/validate -name '*.cc' -exec basename {} \; | sed 's/.cc//g' | sed 's/^/test\/validate\//g')

# Build and run all tests
test: unit perf validate

# Build some examples
example-set: examples/tutorial-1 examples/tutorial-2 examples/hmatComplexityTest examples/hmatWrapperTest examples/ilasParQRSVAAATest \
	examples/inspectMesh examples/kmeansPerfTest examples/nearFieldQRAAA examples/perfTestSVDs examples/printSphere \
	examples/rankRevealDecompTest examples/solveTest examples/testBIEformulations examples/testFunc examples/testNeumann \
	examples/timeComputeKWithRecycling examples/timingsParallel examples/uhmatComplexityTest examples/unfcTestVSC examples/testTACA

# Build some experiments
experiment-set: experiments/benchmarkConstruction experiments/benchmarkMatVec experiments/buildComplexityTest experiments/perfTestACA experiments/plotGNU \
	experiments/plotMatVec experiments/profileMatVec experiments/sphereSolveTest experiments/testSphereConstruction \
	experiments/timingsACA experiments/UNF-meshTest experiments/unfCompression experiments/unfcSphereAnalysis

### INITIALIZATION RECIPES ###
# Create build and bin directories
init:
	@mkdir -p $(DIR_build)/test/unit
	@mkdir -p $(DIR_bin)/test/unit
	@mkdir -p $(DIR_build)/test/perf
	@mkdir -p $(DIR_bin)/test/perf
	@mkdir -p $(DIR_build)/test/validate
	@mkdir -p $(DIR_bin)/test/validate
	@mkdir -p $(DIR_build)/examples
	@mkdir -p $(DIR_bin)/examples
	@mkdir -p $(DIR_build)/experiments
	@mkdir -p $(DIR_bin)/experiments
	@printf "%b" "$(COLOR_BLUE)build and bin directories created.$(NO_COLOR)\n";

### INDIVIDUAL RECIPES ###
# General examples recipe
examples/%: $(DIR_build)examples/%.o
	@mkdir -p $(DIR_bin)$(dir $@)
	@printf "%b" "$(COLOR_BLUE)$(LINK_STRING)\t$(COLOR_CYAN)$(notdir $<)$(NO_COLOR)\n";
	@$(CCOMP) -o $(DIR_bin)$@ $(CFLAGS) $^ $(LFLAGS)

# General experiments recipe
experiments/%: $(DIR_build)experiments/%.o
	@mkdir -p $(DIR_bin)$(dir $@)
	@printf "%b" "$(COLOR_BLUE)$(LINK_STRING)\t$(COLOR_CYAN)$(notdir $<)$(NO_COLOR)\n";
	@$(CCOMP) -o $(DIR_bin)$@ $(CFLAGS) $^ $(LFLAGS)

# General build recipe
$(DIR_build)%.o: $(DIR_scripts)%.cc
	@mkdir -p $(dir $@)
	@printf "%b" "$(COLOR_BLUE)$(COMP_STRING)\t$(COLOR_CYAN)$(notdir $<)$(NO_COLOR)\n";
	@$(CCOMP) -c -o $@ $(CFLAGS) $(IFLAGS) $<


# examples/tutorial-1: $(DIR_build)examples/tutorial-1.o
# 	$(CCOMP) -o $(DIR_bin)$@ $(CFLAGS) $^ $(LFLAGS)

# $(DIR_build)examples/tutorial-1.o: $(DIR_scripts)examples/tutorial-1.cc
# 	$(CCOMP) -c -o $@ $(CFLAGS) $(IFLAGS) $<


# examples/tutorial-2: $(DIR_build)examples/tutorial-2.o
# 	$(CCOMP) -o $(DIR_bin)$@ $(CFLAGS) $^ $(LFLAGS)

# $(DIR_build)examples/tutorial-2.o: $(DIR_scripts)examples/tutorial-2.cc
# 	$(CCOMP) -c -o $@ $(CFLAGS) $(IFLAGS) $<


### TEST RECIPES ###
# Unit tests recipe
$(DIR_run)test/unit/%: $(DIR_bin)test/unit/% FORCE
	@printf "%b" "$(COLOR_BLUE)$(RUN_STRING)\t$(COLOR_CYAN)$(notdir $<)$(NO_COLOR)\n";
	@$<

$(DIR_bin)test/unit/%: $(DIR_build)test/unit/%.o
	@printf "%b" "$(COLOR_BLUE)$(LINK_STRING)\t$(COLOR_CYAN)$(notdir $<)$(NO_COLOR)\n";
	@$(CCOMP) -o $@ $(CFLAGS_unit) $^ $(LFLAGS)

$(DIR_build)test/unit/%.o: $(DIR_scripts)test/unit/%.cc
	@printf "%b" "$(COLOR_BLUE)$(COMP_STRING)\t$(COLOR_CYAN)$(notdir $<)$(NO_COLOR)\n";
	@$(CCOMP) -c -o $@ $(CFLAGS_unit) $(IFLAGS) $<

# Performance tests recipe
$(DIR_run)test/perf/%: $(DIR_bin)test/perf/% FORCE
	@printf "%b" "$(COLOR_BLUE)$(RUN_STRING)\t$(COLOR_CYAN)$(notdir $<)$(NO_COLOR)\n";
	@$<

$(DIR_bin)test/perf/%: $(DIR_build)test/perf/%.o
	@printf "%b" "$(COLOR_BLUE)$(LINK_STRING)\t$(COLOR_CYAN)$(notdir $<)$(NO_COLOR)\n";
	@$(CCOMP) -o $@ $(CFLAGS_perf) $^ $(LFLAGS)

$(DIR_build)test/perf/%.o: $(DIR_scripts)test/perf/%.cc
	@printf "%b" "$(COLOR_BLUE)$(COMP_STRING)\t$(COLOR_CYAN)$(notdir $<)$(NO_COLOR)\n";
	@$(CCOMP) -c -o $@ $(CFLAGS_perf) $(IFLAGS) $<

# Validation tests recipe
$(DIR_run)test/validate/%: $(DIR_bin)test/validate/% FORCE
	@printf "%b" "$(COLOR_BLUE)$(RUN_STRING)\t$(COLOR_CYAN)$(notdir $<)$(NO_COLOR)\n";
	@$<

$(DIR_bin)test/validate/%: $(DIR_build)test/validate/%.o
	@printf "%b" "$(COLOR_BLUE)$(LINK_STRING)\t$(COLOR_CYAN)$(notdir $<)$(NO_COLOR)\n";
	@$(CCOMP) -o $@ $(CFLAGS_validate) $^ $(LFLAGS)

$(DIR_build)test/validate/%.o: $(DIR_scripts)test/validate/%.cc
	@printf "%b" "$(COLOR_BLUE)$(COMP_STRING)\t$(COLOR_CYAN)$(notdir $<)$(NO_COLOR)\n";
	@$(CCOMP) -c -o $@ $(CFLAGS_validate) $(IFLAGS) -MMD $<


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
