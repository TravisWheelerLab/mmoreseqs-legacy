#################################################################################
#		  NAME: MMORE MAKEFILE
# 		AUTHOR: David Rich 
#
#		SOURCE: https://spin.atomicobject.com/2016/08/26/makefile-c-projects/
#################################################################################

# COMPILERS TOOLS:
SHELL 	:= /bin/sh
CC 		:= gcc
CXX 		:= g++
CP 		:= cp 
# COMMANDLINE COMMANDS 
CD 		:= cd 
MKDIR 	:= mkdir
MKDIR_P 	:= $(MKDIR) -p
RM 		:= rm
RMDIR 	:= rmdir

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

# INPUT DIRECTORIES:
SRC_DIR 			?= src/
SCRIPT_DIR 		?= scripts/
TOOL_DIR 		?= tools/
# INPUT SUBDIRECTORIES
SRC_DIRS 		:= $(shell find $(SRC_DIR) -type d)
SCRIPT_DIRS 	:= $(shell find $(SCRIPT_DIR) -type d)

# TOOL DIRECTORIES:
MMSEQS_DIR 		:= $(TOOL_DIR)/mmseqs/
MMSEQS_BIN 		:= $(MMSEQS_DIR)/bin/
HMMER_DIR 		:= $(TOOL_DIR)/hmmer/
HMMER_BIN 		:= $(HMMER_DIR)/bin/

# INPUT FILES:
SRCS 				:= $(shell find $(SRC_DIR) -name *.cpp -or -name *.c -or -name *.s)
HEADERS 			:= $(shell find $(SRC_DIR) -name *.h)
SCRIPTS 			:= $(shell find $(SCRIPT_DIR) -name *.py -or -name *.sh)
TOOLS 			:= $(shell find $(TOOL_DIR) -name hmmbuild -or -name hmmsearch -or -name -mmseqs)

# OUTPUT DIRECTORIES:
BUILD_DIR 		?= build/
INSTALL_DIR 	?= bin/
USER_INSTALL_DIR ?= /usr/bin/
MMORE_BIN 		?= $(BUILD_DIR)
# OUTPUT SUBDIRECTORIES
MMORE_BIN 		?= $(BUILD_DIR)
BUILD_DIRS 		:= $(SRC_DIRS:%=$(BUILD_DIR)/%)

# OUTPUT FILES:
BUILD_OBJS 			:= $(SRCS:%=$(BUILD_DIR)/%.o)
BUILD_DEPS 			:= $(SRCS:%=$(BUILD_DIR)/%.o.d)
BUILD_SCRIPTS 		:= $(SCRIPTS:%=$(BUILD_DIR)/%)
BUILD_TOOLS 		:= $(TOOLS:%=$(BUILD_DIR)/%)
BUILD_PROGRAMS 	:= $(PROGRAMS:%=$(BUILD_DIR)/%)
INSTALL_PROGRAMS	:= $(PROGRAMS:%=$(INSTALL_DIR)/%)

# INCLUDED DIRECTORIES:
INC_DIRS 		:= $(shell find $(SRC_DIRS) -type d)
INC_FLAGS 		:= $(addprefix -I,$(INC_DIRS))
INC_LIBS 		:= $(addprefix -I,$(EASEL_DIR))
# Macro Variables:
INC_DEFS 		:= PROJECT_LOC=$(ROOT_DIR) 			\
						SCRIPT_LOC=$(SCRIPT_DIR) 			\
						MMSEQS_BIN_LOC=$(MMSEQS_BIN) 		\
						HMMER_BIN_LOC=$(HMMER_BIN) 		\
						MMORE_BIN_LOC=$(MMORE_BIN) 		
# INC_DEFS 		:= $(PROGRAMS:%=$(ROOT_DIR)/$(BUILD_DIR)/%)
INC_DEFS 		:= $(addprefix -D,$(INC_DEFS))
# $(info INC_DEFS $(INC_DEFS))

# LINKED LIBRARIES:
CLIBS 			:= -lm -L${EASEL_DIR} -leasel
CFLAGS 			:= $(INC_FLAGS) $(INC_LIBS) -std=gnu99 $(INC_DEFS)
CXXFLAGS 		:= $(INC_FLAGS) $(INC_LIBS) -MMD -MP $(INC_DEFS)

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
$(info *** BUILD-TYPE: $(TYPE) *** )

# PROGRAMS
PROGRAMS = 				\
	$(MAIN_EXEC) 		\
	$(TARGET_EXEC) 	
# OUTPUT PROGRAMS
BUILD_PROGRAMS 	:= $(PROGRAMS:%=$(BUILD_DIR)/%)
INSTALL_PROGRAMS	:= $(PROGRAMS:%=$(INSTALL_DIR)/%)
USER_INSTALL_PROGRAMS := $(PROGRAMS:%=$(USER_INSTALL_DIR)/%)
$(info USER_INSTALL: $(USER_INSTALL_PROGRAMS))

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

# RECIPES
.PHONY: all install user-install submodules submodules-load easel mmseqs hmmer uninstall clean

# Build
all: $(BUILD_PROGRAMS) 
	@echo "*** MAKE ALL ***"

# Program install local file
$(INSTALL_PROGRAMS): $(BUILD_PROGRAMS)
	$(info INSTALL $(@))
	@$(MKDIR_P) $(INSTALL_DIR)
	@$(CP) $(@:$(INSTALL_DIR)/%=$(BUILD_DIR)/%) $(@)

# # Program install local file
# $(USER_INSTALL_PROGRAMS): $(BUILD_OBJS)
# 	$(info USER_INSTALL $(@))
# 	@$(CP) $(@:$(USER_INSTALL_DIR)/%=$(BUILD_DIR)/%) $(@)

# Program build
$(BUILD_PROGRAMS): $(BUILD_SCRIPTS) $(BUILD_OBJS) 
	$(info PROGRAM $@)
	@$(MKDIR_P) $(@D)
	@$(CC) $(BUILD_OBJS) -o $(@)  $(CFLAGS) $(LDFLAGS) $(CLIBS)

# C source
$(BUILD_DIR)/%.c.o: %.c
	$(info CC $@)
	@$(MKDIR_P) $(@D)
	@$(CC) $(CFLAGS) -c $< -o $(@) $(CLIBS) 

# C++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	$(info CXX $@)
	@$(MKDIR_P) $(@D)
	@$(CXX) $(CXXFLAGS) -c $< -o $(@)

# Scripts
$(BUILD_SCRIPTS): $(SCRIPTS)
	$(info SCRIPT $(@))
	@$(MKDIR_P) $(@D)
	@$(CP) $(@:$(BUILD_DIR)/%=%) $(@)

# Tools 
$(BUILD_TOOLS): $(TOOLS)
	$(info TOOLS $(@))
	@$(MKDIR_P) $(@D)
	@$(CP) $(@:$(BUILD_DIR)/%=%) $(@)

# install to binary file
install: $(INSTALL_PROGRAMS)
	@echo "*** MAKE INSTALL ***"

# uninstall program
uninstall:
	@echo "*** MAKE UNINSTALL ***"
	$(info UNINSTALL $(INSTALL_PROGRAMS))
	$(RM) $(INSTALL_PROGRAMS)

# TODO:
# user-install: $(USER_INSTALL_PROGRAMS)
# 	@echo "*** MAKE USER-INSTALL ***"

# user-uninstall: 
# 	@echo "*** MAKE USER-UNINSTALL ***"
# 	$(info UNINSTALL $(USER_INSTALL_PROGRAMS))
# 	$(RM) $(USER_GLOBAL_PROGRAMS)

# install submodules and make them
submodules:
	@echo '*** MAKE SUBMODULE ***'
	make submodules-load
	make easel 
	make mmseqs
	make hmmer

# load all submodules
submodules-load:
	@echo '*** MAKE SUBMODULE-LOAD ***'
	git submodule update --init --recursive

# install easel
easel:
	@echo "*** MAKE EASEL ***"
	git submodule update --init $(ROOT_DIR)/$(EASEL_DIR)
	$(CD) $(ROOT_DIR)/$(EASEL_DIR);  	\
	autoreconf; 								\
	./configure;
	make -C $(ROOT_DIR)/$(EASEL_DIR) --quiet;
	$(CD) $(ROOT_DIR)

# install mmseqs 
mmseqs:
	@echo "*** MAKE MMSEQS ***"
	git submodule update --init $(ROOT_DIR)/$(MMSEQS_DIR)
	$(CD) $(ROOT_DIR);
	$(MKDIR_P) $(ROOT_DIR)/$(MMSEQS_BIN);
	$(CD) $(MMSEQS_DIR);	 														\
	cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. .; 	\
	make -C $(ROOT_DIR)/$(MMSEQS_DIR) --quiet; 							\
	make install -C $(ROOT_DIR)/$(MMSEQS_DIR) --quiet; 					
	$(CD) $(ROOT_DIR)

# install hmmer 
hmmer:
	@echo "*** MAKE HMMER ***"
	git submodule update --init $(ROOT_DIR)/$(HMMER_DIR)
	$(CD) $(ROOT_DIR)
	$(MKDIR) $(HMMER_BIN)
	$(CD) $(HMMER_DIR)
	./configure --prefix=$(ROOT_DIR)/$(HMMER_BIN)
	make -C $(ROOT_DIR)/$(HMMER_DIR) --quiet
	make install -C $(ROOT_DIR)/$(HMMER_DIR) --quiet
	$(CD) $(ROOT_DIR)

# remove build
clean: 
	@echo "*** MAKE CLEAN ***"
	$(info CLEAN $(BUILD_PROGRAMS))
	rm -r $(BUILD_DIR)
	# @$(RM) make_out.txt gmon.out
	# @$(RM) $(BUILD_PROGRAMS)
	# @$(RM) $(BUILD_OBJS) 
	# @$(RM) $(BUILD_DEPS)
	# @$(RM) $(BUILD_SCRIPTS)
	# @$(RMDIR) $(BUILD_DIRS)
	# @$(RMDIR) $(BUILD_DIR)

