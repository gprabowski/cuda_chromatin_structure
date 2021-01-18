###########################################################

## USER SPECIFIC DIRECTORIES ##

# CUDA directory:
CUDA_ROOT_DIR=/usr/local/cuda-11

##########################################################

## CC COMPILER OPTIONS ##

CC=clang++
CC_FLAGS=-std=c++17
CC_LIBS=

##########################################################

## NVCC COMPILER OPTIONS ##

NVCC=nvcc
# NVCC_FLAGS= -arch=sm_75 --std c++17
NVCC_FLAGS= -arch=sm_80 --std c++17
NVCC_LIBS=

# CUDA library directory:
CUDA_LIB_DIR= -L$(CUDA_ROOT_DIR)/lib64
# CUDA include directory:
CUDA_INC_DIR= -I$(CUDA_ROOT_DIR)/include
# CUDA linking libraries:
CUDA_LINK_LIBS= -lcudart

##########################################################

## Project file structure ##

# Source file directory:
SRC_DIR = src

# Object file directory:
OBJ_DIR = bin

# Include header file diretory:
INC_DIR = include

LIB_DIR = src/lib

##########################################################

## Make variables ##

# Target executable name:
EXE = 3dnome_gpu

# Object files:
OBJS = $(OBJ_DIR)/main.o $(OBJ_DIR)/Anchor.o $(OBJ_DIR)/BedRegion.o $(OBJ_DIR)/BedRegions.o $(OBJ_DIR)/Chromosome.o $(OBJ_DIR)/ChromosomesSet.o $(OBJ_DIR)/Cluster.o $(OBJ_DIR)/Density.o $(OBJ_DIR)/Heatmap.o $(OBJ_DIR)/HierarchicalChromosome.o $(OBJ_DIR)/InteractionArc.o $(OBJ_DIR)/InteractionArcs.o $(OBJ_DIR)/LooperSolver.o $(OBJ_DIR)/Settings.o $(OBJ_DIR)/common.o $(OBJ_DIR)/ini.o $(OBJ_DIR)/INIReader.o $(OBJ_DIR)/mtxlib.o $(OBJ_DIR)/rmsd.o

##########################################################

## Compile ##

# Link c++ and CUDA compiled object files to target executable:
$(EXE) : $(OBJS)
	$(CC) $(CC_FLAGS) $(OBJS) -o $@ $(CUDA_INC_DIR) $(CUDA_LIB_DIR) $(CUDA_LINK_LIBS)

# Compile main .cpp file to object files:
$(OBJ_DIR)/%.o : %.cpp
	$(CC) $(CC_FLAGS) -c $< -o $@

# Compile C++ source files to object files:
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp include/%.h
	$(CC) $(CC_FLAGS) -c $< -o $@

#Compile library files
$(OBJ_DIR)/%.o : $(LIB_DIR)/%.cpp $(LIB_DIR)/%.h
	$(CC) $(CC_FLAGS) -c $< -o $@

#Compile library files
$(OBJ_DIR)/%.o : $(LIB_DIR)/%.c $(LIB_DIR)/%.h
	$(CC) $(CC_FLAGS) -c $< -o $@

# Compile CUDA source files to object files:
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cu $(INC_DIR)/%.h
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ $(NVCC_LIBS)

# Clean objects in object directory.
clean:
	$(RM) bin/* *.o $(EXE)





# ORIGINAL MAKEFILE

# SRC=$(wildcard src/*.cpp) $(wildcard src/lib/*.c*)
# OBJ=$(SRC:.cpp=.o)

# CFLAGS=-Wno-write-strings 
# CFLAGSAPP=-l 3dnome -I"./src"  -L"./" -Wl,-rpath="./"

# .PHONY : all

# all: 3dnome 3dnome-app

# 3dnome:
# 		$(CC) $(CFLAGS) -shared $(SRC) -o lib3dnome.so

# 3dnome-app:
# 		$(CC) $(CFLAGS) main.cpp -o 3dnome $(CFLAGSAPP)

# clean:
# 		rm -fr *.o
# 		rm lib3dnome.so
# 		rm 3dnome
