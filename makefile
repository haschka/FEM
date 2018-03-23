# Hexahedron mesh finite element solver. 
#
# compilier and compilier flags

#CC = x86_64-w64-mingw32-gcc
CC = gcc

#CFLAGS = -g -mssse3
#CFLAGS = -g -mavx2 -mfma
#CFLAGS = -O3 -mssse3 -march=core2 -flto -fopt-info-vec-optimized -fomit-frame-pointer
#CFLAGS = -O3 -mssse3 -march=native -flto -fopt-info-vec-optimized -fprofile-arcs -ftest-coverage -g 
CFLAGS = -O3 -march=native -flto -fopt-info-vec-optimized -fomit-frame-pointer 
#CFLAGS = -O3 -mtune=750 -mcpu=750 -flto -fomit-frame-pointer

OUT = hex_fem
DEBUG = 

MATH= -lm

# Comment out on OSX math is part of the lsystem

all: mesh_from_file \
 shape_func \
 matrix_math \
 node_shuffle \
 sparse_matrix \
 stiffness_matrix_for_element \
 global_stiffness_matrix \
 post_processing

	$(CC) $(CFLAGS) $(DEBUG) $(LAPACK) $(CHOLMOD) $(AMD) $(MATH) \
	-o $(OUT) hex_fem_solver.c \
					mesh_from_file.o \
					shape_func.o \
					matrix_math.o \
					node_shuffle.o \
					sparse_matrix.o \
					stiffness_matrix_for_element.o \
					global_stiffness_matrix.o \
					post_processing.o

mesh_from_file:
	$(CC) $(CFLAGS) $(DEBUG) -c mesh_from_file.c

shape_func:
	$(CC) $(CFLAGS) $(DEBUG) -c shape_func.c

matrix_math:
	$(CC) $(CFLAGS) $(DEBUG) -c matrix_math.c

node_shuffle:
	$(CC) $(CFLAGS) $(DEBUG) -c node_shuffle.c

sparse_matrix:
	$(CC) $(CFLAGS) $(DEBUG) -c sparse_matrix.c

stiffness_matrix_for_element: shape_func matrix_math mesh_from_file
	$(CC) $(CFLAGS) $(DEBUG) -c stiffness_matrix_for_element.c

global_stiffness_matrix: stiffness_matrix_for_element
	$(CC) $(CFLAGS) $(DEBUG) -c global_stiffness_matrix.c

post_processing: global_stiffness_matrix
	$(CC) $(CFLAGS) $(DEBUG) -c post_processing.c
