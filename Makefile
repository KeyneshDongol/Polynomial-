###############################################################################
#                                                                             #
#                  Makefile for projects utilising TORTOISE                   #
#                                                                             #
##################################### Usage ###################################
#                                                                             #
# 1) Place this makefile in the same folder as the main.cpp file              #
# 2) Specify the location of the dependency files                             #
# 3) Edit executable name, compiler flags and the like (optional)             #
# 4) Run make all                                                             #
#                                                                             #
# PS: Do not forget to specify the path to GNUPlot in the Configuration file  #
#                                                                             #
############################## Editable settings ##############################

# Location of the library Eigen
EIGE_DIR := /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3
# Location of the project main file and subdirectories
USER_DIR := .
# Name of the main file
EXE      := main
# Compiler
CXX      := /opt/homebrew/Cellar/llvm/16.0.6/bin/clang++
# Compilation flags
FLAGS    := -std=c++20 -O3 -flto -fopenmp
#If you want to debug, uncomment the following line:
#FLAGS    := -glldb -g3 -std=c++20 -O0 -fopenmp


###############################################################################

####################### Do not edit below this point ##########################

# sources, objects and headers from the base directory
SRCS_USR := main.cpp
OBJS := main.o
DEPS := main.d

# folders to include with -I
INC_FLAGS := -I Polynomial -I $(EIGE_DIR)

# c++ flags, -I and automatic dependency generator
CPPFLAGS := $(FLAGS) $(INC_FLAGS) -MMD -MP

########################### Makefile rules start here ########################

default :
	@echo 'Typing "make" with no arguments displays this message.'
	@echo 'Specific targets for compiling are:'
	@echo '  all       compiles the whole project'
	@echo '  clean     removes executable, objects, archive and dependencies'

.PHONY: clean

all : $(EXE)

$(EXE) : $(OBJS)
	$(CXX)  $(OBJS) -o $@ $(LDFLAGS) -L/opt/homebrew/Cellar/llvm/16.0.6/lib/ -L/usr/local/include -fopenmp -lomp

-include $(DEPS)

clean :
	@rm -f $(OBJS) $(EXE) $(DEPS)

