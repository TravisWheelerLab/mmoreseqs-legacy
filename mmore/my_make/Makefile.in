#################################################################################
#		  NAME: TEMPLATE MAKEFILE
# 		AUTHOR: David Rich 
#
#		SOURCE: https://spin.atomicobject.com/2016/08/26/makefile-c-projects/
#################################################################################

# COMPILER TOOLS:
SHELL 	= /bin/sh
CC 		= gcc
CXX 		= g++
# OTHER TOOLS:
CP 		= cp

# ROOT FOLDER:
ROOT_DIR ?= $(shell pwd)
# $(info **** ROOT_DIR = $(ROOT_DIR) **** )

# PACKAGE/BUILD INFORMATION:
PACKAGE 				= FB-PRUNER 
PACKAGE_VERSION 	= 0.1A
PACKAGE_TARNAME 	= mmore

# FINAL EXECUTABLE:
TARGET_EXEC 	?= mmore
MAIN_EXEC 		?= mmore

# LIBRARY DIRECTORIES: 
LIB_DIR 			?= ./lib
EASEL_DIR 		?= $(LIB_DIR)/easel

# main tool directory (need this for scripts)
MMORE_BIN 		?= build/

# TOOL DIRECTORIES:
TOOL_DIR 		?= tools/
MMSEQS_DIR 		?= $(TOOL_DIR)/MMseqs2
MMSEQS_BIN 		?= $(MMSEQS_DIR)/bin
HMMER_DIR 		?= $(TOOL_DIR)/hmmer
HMMER_BIN 		?= $(HMMER_DIR)/bin

# INPUT/OUTPUT DIRECTORIES:
INSTALL_DIR 	?= bin/
BUILD_DIR 		?= build/
SRC_DIRS 		?= src/
SCRIPT_DIRS 	?= scripts/

# SOURCES:
# get all source code from source directory
SRCS 				:= $(shell find $(SRC_DIRS) -name *.cpp -or -name *.c -or -name *.s)
# get all scripts from script directory
SCRIPTS 			:= $(shell find $())
# for each source file, create an object (.o)
OBJS 				:= $(SRCS:%=$(BUILD_DIR)/%.o)
# for each object file create a dependency file (.d)
DEPS 				:= $(OBJS:.o=.d)

# INCLUDED DIRECTORIES:
INC_DIRS 		:= $(shell find $(SRC_DIRS) -type d)
# Include flags:
INC_FLAGS 		:= $(addprefix -I,$(INC_DIRS))

# Include library archives
INC_LIBS 		:= -I$(EASEL_DIR)
# Macro Variables:
INC_VARS 		:= -DPROJECT_LOC=$(ROOT_DIR) 			\
						-DMMSEQS_BIN_LOC=$(MMSEQS_BIN) 	\
						-DHMMER_BIN_LOC=$(HMMER_BIN) 		\
						-DMMORE_BIN_LOC=$(MMORE_BIN) 

# LINKED LIBRARIES:
# suggested flags: -lm -L(easel library location) -leasel
CLIBS 		:= -lm -L${EASEL_DIR} -leasel
# standard C compiler flags to be run (non-build-specific)
# suggested flags: -std=gnu99 -g -DROOT_DIR=$(ROOT_DIR)
CFLAGS 		:= $(INC_FLAGS) $(INC_LIBS) -std=gnu99 $(INC_VARS)
# standard C pre-processor compiler flags to be run (project currently contains no current c++ files)
# suggested flags: -MMD -MP -O3
CXXFLAGS 	:= $(INC_FLAGS) $(INC_LIBS) -MMD -MP $(INC_VARS)

# BUILD-SPECIFIC FLAGS: 
# -g (DEBUGGING INFO FOR VALGRIND)
# -pg (PROFILING INFO FOR GPROF)
# suggested flags: -g (line numbers) -pg (gprof) -DDEBUG=1 -fsanitize=address
DEBUGFLAGS 		:= -g -pg -fsanitize=address -DDEBUG=1 
# flags for alpha version
ALPHAFLAGS 		:= -g -O1 -fsanitize=address
# flags for alpha version
VALFLAGS 		:= -g -O0 
# flags for release version
RELEASEFLAGS 	:= -O3

# SELECT BUILD-SPECIFIC FLAGS: 
# special flags for DEBUG
ifeq ($(BUILD), DEBUG)
	TYPE			=  DEBUG
	CFLAGS		+= $(DEBUGFLAGS)
	CPPFLAGS		+=	$(DEBUGFLAGS)
	TARGET_EXEC	=	mmore-DEBUG
endif
# special flags for ALPHA 
ifeq ($(BUILD), ALPHA)
	TYPE			=	ALPHA
	CFLAGS		+=	$(ALPHAFLAGS)
	CPPFLAGS		+=	$(ALPHAFLAGS)
	TARGET_EXEC	=	mmore-ALPHA
endif
# special flags for ALPHA 
ifeq ($(BUILD), VALGRIND)
	TYPE			=	VALGRIND
	CFLAGS		+=	$(VALFLAGS)
	CPPFLAGS		+=	$(VALFLAGS)
	TARGET_EXEC	= 	mmore-VALGRIND
endif
# special flags for RELEASE 
BUILD ?= RELEASE
ifeq ($(BUILD), RELEASE)
	TYPE			=	RELEASE
	CFLAGS		+=	$(RELEASEFLAGS)
	CPPFLAGS		+=	$(RELEASEFLAGS)
	TARGET_EXEC	=	mmore-RELEASE
endif
# build type
$(info *** BUILD-TYPE: ${TYPE} *** )

# PROGRAMS
PROGRAMS = \
	${BUILD_DIR}/${TARGET_EXEC} \

# $(info *** PROGRAMS: ${PROGRAMS} *** )
# $(info *** CFLAGS: ${CFLAGS} *** )

# SILENCE UNNECCESSARY OUTPUT
ECHO_SUBDIR0 = +${MAKE} -C #space separator after -c
ECHO_SUBDIR1 = 
ifndef V
	QUIET_CC      = echo 'CC $@'; ${CC}
	QUIET_CXX     = echo 'CXX $@'; ${CXX}
	QUIET_GEN     = echo 'GEN $@'; ${GEN}
	QUIET_AR      = echo 'AR $@'; ${AR}
	QUIET_SUBDIR0 = +@subdir=
	QUIET_SUBDIR1 = ; echo 'SUBDIR $$subdir' \
		        		echo '${MAKE} -s -C $$subdir'
endif

.PHONY: all force install clean

# full install
all: ${PROGRAMS}

# C executable
${PROGRAMS}: $(OBJS)
	$(info CC $@)
	@$(CC) $(OBJS) -o $@  $(CFLAGS) $(LDFLAGS) $(CLIBS)

# C source
$(BUILD_DIR)/%.c.o: %.c
	$(info CC $@)
	@$(MKDIR_P) $(dir $@)
	@$(CC) $(CFLAGS) -c $< -o $@ $(CLIBS) 

# C++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	$(info CXX $@)
	@$(MKDIR_P) $(dir $@)
	@$(CXX) $(CXXFLAGS) -c $< -o $@

# force remake
force:
	echo '*** MAKE FORCE-ALL ***'

# install submodules and make them
submodules:
	echo '*** MAKE SUBMODULE ***'
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