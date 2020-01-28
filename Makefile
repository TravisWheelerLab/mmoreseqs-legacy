#################################################################################
#
#		TEMPLATE MAKEFILE
#
#		SOURCE: https://spin.atomicobject.com/2016/08/26/makefile-c-projects/
#
#################################################################################

# define C compiler to use
CC = gcc-8

# final executable
TARGET_EXEC ?= cloud_fwdbck.exe

# input/output directories
BUILD_DIR ?= ./build
SRC_DIRS ?= ./src

# get all source code from source directory
SRCS := $(shell find $(SRC_DIRS) -name *.cpp -or -name *.c -or -name *.s)
# for each source file, create an object (.o)
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
# for each object file create a dependency file (.d)
DEPS := $(OBJS:.o=.d)

# find all included dependency files
INC_DIRS := $(shell find $(SRC_DIRS) -type d)
# create a compiler flag for each included dependency
INC_FLAGS := $(addprefix -I,$(INC_DIRS))

# Any C linked libraries I want to include
CLIBS := -lm
# standard C compiler flags to be run
# suggested flags: -Wall -g -O2 -DDEBUG
CFLAGS := $(INC_FLAGS) -g -std=c99
# standard C++ compiler flags to be run
CPPFLAGS ?= $(INC_FLAGS) -MMD -MP

# compile each object file before executable
$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS) $(CLIBS)

# assembly
$(BUILD_DIR)/%.s.o: %.s
	$(MKDIR_P) $(dir $@)
	$(AS) $(ASFLAGS) -c $< -o $@

# c source
$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(CLIBS) -c $< -o $@

# c++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@


.PHONY: clean

clean:
	$(RM) -r $(BUILD_DIR)

-include $(DEPS)

MKDIR_P ?= mkdir -p

#
