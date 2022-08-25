#################################################################################
#  - FILE: MMOREseqs Makefile
#  - DESC: This is root Makefile.  This will run CMake in proper subdirectory.
#################################################################################

# Compilers
SHELL 		:= /bin/sh
CC 				:= gcc
CXX 			:= g++
# COMMANDLINE COMMANDS 
CD 				:= cd 
MKDIR 		:= mkdir
MKDIR_P 	:= $(MKDIR) -p
RM 				:= rm
RMDIR 		:= rmdir

# root directory.
ROOT_DIR ?= $(shell pwd)
# Input directories
SOURCE_DIR 			?= src/
SCRIPT_DIR 			?= scripts/
# Input subdirectories
SOURCE_DIRS 		:= $(shell find $(SOURCE_DIR) -type d)
SCRIPT_DIRS 		:= $(shell find $(SCRIPT_DIR) -type d)
# Input files
SOURCES 				:= $(shell find $(SOURCE_DIR) -name *.cpp -or -name *.c -or -name *.s)
HEADERS 				:= $(shell find $(SOURCE_DIR) -name *.hpp -or -name *.h)
SCRIPTS 				:= $(shell find $(SCRIPT_DIR) -name *.py -or -name *.sh)
TOOLS 					:= $(shell find $(TOOL_DIR) -name hmmbuild -or -name hmmsearch -or -name -mmseqs)

# get commit hash of this build.
BUILD_HASH 			:= $(shell git rev-parse HEAD)

# Build directories
BUILD_RELEASE_DIR 	:= build
BUILD_DEBUG_DIR 		:= build-debug

# Build-Specific Flags
# -g (Debugging symbols)
# -pg (Profiling info)
# flags for release version
RELEASE_C_FLAGS 	:= -O3 -Wall
# suggested flags: -g (line numbers) -pg (gprof) -DDEBUG=1 -fsanitize=address
DEBUG_C_FLAGS 		:= -g -O1 -pg -fsanitize=address
# flags for alpha version
VALGRIND_C_FLAGS 	:= -g -O0 

# RECIPES
.PHONY: all build-release build-debug build-valgrind build-test test format clean clean-release clean-debug

default:
	$(MAKE) build-release

build-cli:
	@echo "*** BUILD CLI ***"
	

build-release:
	@echo "*** BUILD RELEASE ***"
	@mkdir -p $(BUILD_RELEASE_DIR)
	@cmake -S . -B $(BUILD_RELEASE_DIR) -DCMAKE_BUILD_TYPE=RELEASE -DSET_BUILD_HASH=$(BUILD_HASH)
	@cd $(BUILD_RELEASE_DIR) && \
		make CFLAGS="$(RELEASE_C_FLAGS)" CXXFLAGS="$(RELEASE_C_FLAGS)"

build-debug:
	@echo "*** BUILD DEBUG ***"
	@mkdir -p $(BUILD_DEBUG_DIR)
	@cmake -S . -B $(BUILD_DEBUG_DIR) -DCMAKE_BUILD_TYPE=DEBUG -DSET_BUILD_HASH=$(BUILD_HASH)
	@cd $(BUILD_DEBUG_DIR) && \
		make CFLAGS="$(DEBUG_C_FLAGS)" CXXFLAGS="$(DEBUG_C_FLAGS)" 

build-valgrind:
	@echo "*** BUILD VALGRIND ***"
	@mkdir -p $(BUILD_DEBUG_DIR) 
	@cmake -S . -B $(BUILD_DEBUG_DIR) -DCMAKE_BUILD_TYPE=DEBUG -DSET_BUILD_HASH=$(BUILD_HASH)
	@cd $(BUILD_DEBUG_DIR) && \
		make CFLAGS="$(VALGRIND_C_FLAGS)" CXXFLAGS="$(VALGRIND_C_FLAGS)"

build-test:
	@echo '*** BUILD-TEST ***'
	$(MAKE) build-release 
	@echo ''
	$(MAKE) test

test:
	@echo '*** TEST ***'
	./$(BUILD_RELEASE_DIR)/bin/mmoreseqs_test

# retrieve all submodules
get-submodules:
	@echo '*** GET SUBMODULE ***'
	git submodule update --init --recursive 

# code formatting 
format: $(SOURCES)
	@echo '*** FORMAT ***'
	clang-format -i -style=file $(SOURCES) $(HEADERS)

clean: 
	@echo "*** CLEAN ***"
	$(MAKE) clean-release 
	$(MAKE) clean-debug

clean-release: 
	@echo "*** CLEAN-RELEASE ***"
	@cd $(BUILD_RELEASE_DIR) && \
		make clean 

clean-debug:
	@echo "*** CLEAN-DEBUG ***"
	@cd $(BUILD_DEBUG_DIR) && \
		make clean 
