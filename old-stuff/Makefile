#################################################################################
#		  NAME: TEMPLATE MAKEFILE
# 		AUTHOR: David Rich (modified from source)
#
#		SOURCE: https://spin.atomicobject.com/2016/08/26/makefile-c-projects/
#################################################################################

# define C compiler to use
SHELL = /bin/sh
CC = gcc
CXX = g++

# location of root folder of makefile
ROOT_DIR ?= $(shell pwd)
$(info **** ROOT_DIR = $(ROOT_DIR) **** )

# package information 
PACKAGE = FB-PRUNER 
PACKAGE_VERSION = 0.1A
PACKAGE_TARNAME = mmore

# final executable
TARGET_EXEC ?= mmore
MAIN_EXEC 	?= mmore

# library directories 
LIB_DIR 		?= ./lib
EASEL_DIR 	?= $(LIB_DIR)/easel

# main tool directory (need this for scripts)
MMORE_BIN 	?= ./build

# other tools directories
TOOL_DIR 	?= ./tools
MMSEQS_DIR 	?= $(TOOL_DIR)/MMseqs2
MMSEQS_BIN 	?= $(MMSEQS_DIR)/bin
HMMER_DIR 	?= $(TOOL_DIR)/hmmer
HMMER_BIN 	?= $(HMMER_DIR)/bin


# input/output directories
INSTALL_DIR ?= ./bin
BUILD_DIR 	?= ./build
SRC_DIRS 	?= ./src

# get all source code from source directory
SRCS 	:= $(shell find $(SRC_DIRS) -name *.cpp -or -name *.c -or -name *.s)
# for each source file, create an object (.o)
OBJS 	:= $(SRCS:%=$(BUILD_DIR)/%.o)
# for each object file create a dependency file (.d)
DEPS 	:= $(OBJS:.o=.d)

# find all included dependency files
INC_DIRS 	:= $(shell find $(SRC_DIRS) -type d)
# create a compiler flag for each included dependency
INC_FLAGS 	:= $(addprefix -I,$(INC_DIRS))

# include library archives
INC_LIBS := -I$(EASEL_DIR)

# variables declarations
INC_VARS := -DPROJECT_LOC=$(ROOT_DIR) -DMMSEQS_BIN_LOC=$(MMSEQS_BIN) -DHMMER_BIN_LOC=$(HMMER_BIN) -DMMORE_BIN_LOC=$(MMORE_BIN) 

# Any C linked libraries I want to include
# suggested flags: -lm -L(easel library location) -leasel
CLIBS 			:= -lm -L${EASEL_DIR} -leasel
# standard C compiler flags to be run (non-build-specific)
# suggested flags: -std=gnu99 -g -DROOT_DIR=$(ROOT_DIR)
CFLAGS 			:= $(INC_FLAGS) $(INC_LIBS) -std=gnu99 $(INC_VARS)
# standard C pre-processor compiler flags to be run (project currently contains no current c++ files)
# suggested flags: -MMD -MP -O3
CPPFLAGS 		?= $(INC_FLAGS) $(INC_LIBS) -MMD -MP $(INC_VARS)
# flags for debugging and profiling
# suggested flags: -g -pg -DDEBUG=1 -fsanitize=address
DEBUGFLAGS 		:= $(INC_FLAGS) $(INC_LIBS) -g -pg -DDEBUG=1 -fsanitize=address
# flags for alpha version
ALPHAFLAGS 		:= $(INC_FLAGS) $(INC_LIBS) -fsanitize=address
# flags for release version
RELEASEFLAGS 	:= $(INC_FLAGS) $(INC_LIBS) -O3

$(info ************  BUILD FB_PRUNER  ************)
# special flags for DEBUG
ifeq ($(BUILD), DEBUG)
	CFLAGS+=$(DEBUGFLAGS)
	CPPFLAGS+=$(DEBUGFLAGS)
	TARGET_EXEC=mmore-DEBUG
endif
# special flags for ALPHA 
ifeq ($(BUILD), ALPHA)
# 	$(info RELEASE FLAGGED)
	CFLAGS+=$(ALPHAFLAGS)
	CPPFLAGS+=$(RELEASEFLAGS)
	TARGET_EXEC=mmore-ALPHA
endif
# special flags for RELEASE 
ifeq ($(BUILD), RELEASE)
# 	$(info RELEASE FLAGGED)
	CFLAGS+=$(RELEASEFLAGS)
	CPPFLAGS+=$(RELEASEFLAGS)
	TARGET_EXEC=mmore-RELEASE
endif

# compile each object file before executable
$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(CC) $(OBJS) -o $@  $(CFLAGS) $(LDFLAGS) $(CLIBS)

# assembly
$(BUILD_DIR)/%.s.o: %.s
	$(MKDIR_P) $(dir $@)
	$(AS) $(ASFLAGS) -c $< -o $@

# c source
$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@ $(CLIBS) 

# c++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY: all install clean

# full install
all:
	# compile each object file before executable
	$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
		$(CC) $(OBJS) -o $@  $(CFLAGS) $(LDFLAGS) $(CLIBS)

	# assembly
	$(BUILD_DIR)/%.s.o: %.s
		$(MKDIR_P) $(dir $@)
		$(AS) $(ASFLAGS) -c $< -o $@

	# c source
	$(BUILD_DIR)/%.c.o: %.c
		$(MKDIR_P) $(dir $@)
		$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@ $(CLIBS) 

	# c++ source
	$(BUILD_DIR)/%.cpp.o: %.cpp
		$(MKDIR_P) $(dir $@)
		$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# pull submodules and make them
submodules:
	git submodule update --recursive
	make easel 
	make mmseqs

# install easel
easel:
	cd $(ROOT_DIR)
	cd $(EASEL_DIR)
	autoconf
	./configure
	$(MAKE) -C $(EASEL_DIR)

# install mmseqs 
mmseqs:
	cd $(ROOT_DIR)
	cd $(MMSEQS_DIR)
	cmake

# remove build
clean:
	$(RM) -r $(BUILD_DIR) make_out.txt gmon.out

-include $(DEPS)

MKDIR_P ?= mkdir -p

install:
	# mkdir $(INSTALL_DIR)
	cp -f $(BUILD_DIR)/$(TARGET_EXEC) $(BUILD_DIR)/$(MAIN_EXEC)
	# cp -f $(BUILD_DIR)/$(TARGET_EXEC) $(INSTALL_DIR)/

uninstall:
	rm 