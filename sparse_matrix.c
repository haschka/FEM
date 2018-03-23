/*! \file sparse_matrix.c
 *  \brief The implemenation of sparse_matricies.
 *  This file contains all functions necessary to handle the sparse matricies
 *  used during this finite element process. In this case the global
 *  stiffness matrix in full and reduced form. The main feature is the 
 *  bisection base search which allows us to use sparse matricies at 
 *  reasonable runtime.
 */

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>
#include"hex_fem_solver.h"

/*! \brief Matrix sorting (ordering) helper structure. 
 *  Auxilary matrix handling structure that contains the global index
 *  for a value going from 0 to range*range-1. and the corresponding
 *  value if it is used. No values and indicies for values = zero are stored
 *  This structure is mainly used for sorting.
 */
typedef struct {
  long long global_index; /*!< global index of 0 to range*range-1 */
  double value; /*!< value at index */
} sparse_matrix_order_handle;

/*! \brief a function that allocates a binary ( index, hint, eye) 
 *         matrix of range range.
 *  This function allocates the memory for a binary index matrix that might
 *  be used for searching values in a sparse matrix efficiently later on.
 *  \param range range of the index matrix.
 */
char* create_binary_index_matrix(long long range) {
  char* matrix = (char*)malloc((range*range)/(sizeof(char)*8)+1);
  memset(matrix,0,(range*range)/(sizeof(char)*8)+1);
  return(matrix);
}

/*! \brief a function that sets a value in a binary index ( hint, eye ) matrix.
 *  This function sets a value in a binary index matrix, i.e. makes sure
 *  that at the global index of the matrix the value is set to 1. 
 *  \param bin_index_matrix The binary index matrix to operate on.
 *  \param index The global index (between 0 and range*range-1)
 *               where to set the value to 1. 
 */
void set_value_in_binary_index_at_index(char* bin_index_matrix,
					long long index) {

  char small_index;
  
  int module = index%(sizeof(char)*8);
  int location = index/(sizeof(char)*8);
  
  small_index = (char)1 << module;

  bin_index_matrix[location] = bin_index_matrix[location] | small_index;
}

/*! \brief A function that retrieves a value from a binary index ( hint, eye )
 *         matrix at a specified index. 
 *  This function retrieves a value from a binary index matrix at the 
 *  specified global index. (returns either 0 or 1 as stored at the index).
 *  \param bin_index_matrix The binary index matrix to operate on.
 *  \param index The global index (between 0 and range*range-1)
 *               from where to retrieve the value.  
 */
int get_value_in_binary_index_at_index(char* bin_index_matrix,
				       long long index) {

  char small_index;
  int module = index%(sizeof(char)*8);
  int location = index/(sizeof(char)*8);

  small_index = (char)1 << module;

  return (small_index == (bin_index_matrix[location] & small_index));
} 

/*! \brief A comparison function in order to sort a sparse matrix using qsort.
 *  This function allows for sorting of the values in a sparse matrix using
 *  the qsort algorith, and sparse matricies that are store in 
 *  sparse_matrix_order_handle form. The function returns 0 if the 
 *  global indicies of p and q are the same, 1 if the global index of p > q,
 *  -1 if the global index of p < q. 
 *  \param p a pointer to a sparse_matrix_order_handle offset (const void here)
 *  \param q a pointer to a sparse_matrix_order_handle offset (const void here)
 */
int sparse_matrix_order_handle_compare(const void *p, const void *q) {
  sparse_matrix_order_handle x = *(const sparse_matrix_order_handle*)p;
  sparse_matrix_order_handle y = *(const sparse_matrix_order_handle*)q;

  int retval;

  if ( x.global_index == y.global_index ) {
    retval = 0;
  }
  if ( x.global_index > y.global_index ) {
    retval = 1;
  }
  if ( x.global_index < y.global_index) {
    retval = -1;
  }
  return (retval);
}

/*! \brief A function that orders sparse matricies by its index.
 *  This function orders a sparse matrix by its index. This is important
 *  for efficient random access value retrieval.
 */
void order_sparse_matrix(sparse_matrix s_mat) {
  long long i;

  sparse_matrix_order_handle* handle =
    (sparse_matrix_order_handle*)malloc(sizeof(sparse_matrix_order_handle)
					*s_mat.non_zeros);

  for(i=0;i<s_mat.non_zeros;i++) {
    handle[i].global_index =
      (long long)s_mat.row_indicies[i]*(long long)s_mat.range
      +(long long)s_mat.col_indicies[i];

    handle[i].value = s_mat.values[i];
  }

  qsort(handle,
	(size_t)s_mat.non_zeros,
	sizeof(sparse_matrix_order_handle),
	sparse_matrix_order_handle_compare);

  for(i=0;i<s_mat.non_zeros;i++) {
    s_mat.row_indicies[i] =
      (long long)handle[i].global_index/(long long)s_mat.range;
    s_mat.col_indicies[i] =
      (long long)handle[i].global_index%(long long)s_mat.range;
    s_mat.values[i] = handle[i].value;
  }
  free(handle);
}

/*! \brief function that gernerates a sparse matrix from an index matrix
 *  This function reserves and prepares the sparse matrix from an index 
 *  matrix, it stores rank and collom but not values that have to be filled 
 *  in later on.
 *  \param index_matrix index matrix that hints where this sparse matrix 
 *                      will have non zero values
 *  \param range the index matrix rank ( rank of the sparse matrix is 3*range )
 *               as the index matrix stores indicies for entire nodes.
 */
sparse_matrix generate_sparse_matrix_from_index_matrix(char* index_matrix,
						       int range) {

  int i,j,k,l,m;
  int i_three, j_three, k_nine;
  long long non_zeros;
  
  sparse_matrix s_mat;
  s_mat.range = range*3;

  non_zeros =
    number_of_non_zero_elements_in_index_matrix(index_matrix, range);

  s_mat.row_indicies = (int*)malloc(sizeof(int)*non_zeros*9);
  s_mat.col_indicies = (int*)malloc(sizeof(int)*non_zeros*9);

  s_mat.values = (double*)malloc(sizeof(double)*non_zeros*9);
  memset(s_mat.values,0,sizeof(double)*non_zeros*9);
  
  k=0;
  for(i=0;i<range;i++) {
    i_three = i*3;
    for(j=0;j<range;j++) {
      if(index_matrix[((long long)i*(long long)range)+(long long)j] == 1) {
	j_three = j*3;
	for(l=0;l<3;l++) {
	  for(m=0;m<3;m++) {
	    s_mat.col_indicies[k] = i_three+l;
	    s_mat.row_indicies[k] = j_three+m;
	    k++;
	  }
	}
      }
    }
  }
#ifdef _debug
  if (k != non_zeros*9) {
    printf("DEBUG! Counters are wrong: k != non_zeros \n");
    printf("DEBUG! @ - generate_sparse_matrix_from_index_matrix() \n");
  }
#endif
  s_mat.has_eye = 0;
  s_mat.non_zeros = non_zeros*9;
  return(s_mat);
}

/*! \brief a function that frees the memory of a sparse matrix 
 *  This function destroys a sparse matrix and frees it memory
 */
void free_sparse_matrix(sparse_matrix s_mat) {
  free(s_mat.row_indicies);
  free(s_mat.col_indicies);
  free(s_mat.values);
  if (s_mat.has_eye == 1) {
    free(s_mat.eye);
  }
}

/*! \brief a function that generates a binray matrix that hints weather 
 *         a value is stored at a certain location in a sparse matrix or not
 *  A function that generates a binray matrix that is 1 where a value is stored
 *  and 0 where no value is stored. This shall speed up finding values in the
 *  matrix and allow for quick 0 returns if not value is stored at a specific
 *  location
 */ 
char* create_sparse_eye(sparse_matrix s_mat) {

  int i;
  char* eye = create_binary_index_matrix((long long)s_mat.range);

  for(i=0;i<s_mat.non_zeros;i++) {
    set_value_in_binary_index_at_index(eye,
				       (long long)s_mat.row_indicies[i]
				       *(long long)s_mat.range
				       +(long long)s_mat.col_indicies[i]);
  }
  return(eye);
}

/*! \brief A function to find the index of a specific value in a sparse matrix.
 *  The function returns the global array index for a value stored at column 
 *  col and at row row. This function implements a bisection algorithm 
 *  and hence requires an ordered ( sorted ) matrix. Bisection speeds up 
 *  the overall program by up to two magnitudes.
 *  \param s_mat the sparse matrix to get the global index from.
 *  \param col the column index. 
 *  \param row the row index
 */
int get_sparse_matrix_index_ordered(sparse_matrix s_mat, int col, int row) {

  int i;

  long long index = (long long)row*(long long)s_mat.range+(long long)col;
  long long current_index;
  int found = 0;
  int bifourcation_index;
  int length = s_mat.non_zeros/2;
  int s_point = length;
  int s_point_new;
  int delta;
  int upper_bound = s_mat.non_zeros-1;

  // go and bifurcate
  
  while(length > 8) {

    delta = length/2+1; 
    current_index =
      (long long)s_mat.row_indicies[s_point]*(long long)s_mat.range+
      (long long)s_mat.col_indicies[s_point];

    // backwards
    if ( current_index > index ) {
      s_point = s_point - delta;
      if ( __builtin_expect(s_point < 0,0)) {
	s_point = 0;
      }
    }
    // forwards
    if ( current_index < index ) {
      s_point = s_point + delta;
      if ( __builtin_expect(s_point > upper_bound,0)  ) {
	s_point = upper_bound;
      }
    }
    // We stumbled on the point
    if ( current_index == index) {
      return (s_point);
    }
    length = delta;
  }

  if ( (long long)s_mat.row_indicies[s_point]
       *(long long)s_mat.range+(long long)s_mat.col_indicies[s_point]
       == index ) {
    return (s_point);
  }
  
  // final search in a small portion
  if ( (long long)s_mat.row_indicies[s_point]
       *(long long)s_mat.range+(long long)s_mat.col_indicies[s_point]
       > index ) {
    // avoid overrun
    if(s_point - length < 0) {
      length = s_point;
    }
    for(i=0;i<length;i++) {
      if ( row == s_mat.row_indicies[s_point-i] &&
	   col == s_mat.col_indicies[s_point-i] ) {
	return(s_point-i);
      }
    }
  } else {
    if (s_point + length > s_mat.non_zeros) {
      length = s_mat.non_zeros - s_point;
    }
    // avoid overrun
    for(i=0;i<length;i++) {
      if ( row == s_mat.row_indicies[s_point+i] &&
	   col == s_mat.col_indicies[s_point+i] ) {
	return(s_point+i);
      }
    }
  }
  // error no index could be found
  printf("Index not found in - @ get_sparse_matrix_index_ordered () \n");
  exit(-1);
}

/*! \brief a function that allows to increment a value in a sparse matrix.
 *  This function allows to incriment a value in a sparse matrix. As sparse 
 *  matrix indexed values are initialized with zero this function is also 
 *  used to set the values of the matrix. 
 *  \param s_mat sparse matrix to operate on.
 *  \param col column index to operate on.
 *  \param row row index to operate on.
 *  \param value value to add to the value stored at the indicies in the matrix.
 */
void add_to_value_in_sparse_matrix(sparse_matrix s_mat,
				   int col,
				   int row,
				   double value) {

  int i = get_sparse_matrix_index_ordered(s_mat,col, row);
  //int i = fast_sparse_matrix_index(s_mat,col,row);
  s_mat.values[i] += value;
}

/*! \brief a function that allows you to obtain a value from a sparse matrix
 *         at column col and row row.
 *  This function allows you to retrieve a value from the sparse matrix by 
 *  column and row index.
 *  \param s_mat sparse matrix to retrieve the value from. 
 *  \param col column of the value to retrieve.
 *  \param row row of the value to retrieve.
 */
double get_value_from_sorted_sparse_matrix(sparse_matrix s_mat,
					   int col, int row) {

  int i = get_sparse_matrix_index_ordered(s_mat,col, row);

#ifdef _debug  
  int j = fast_sparse_matrix_index(s_mat,col,row);

  if ( i != j ) {
    printf("Falal error in finding correct indicies in \n");
    printf("-@ get_value_from_sorted_sparse_matrix\n");
  }
#endif
  return(s_mat.values[i]);
}

/*! \brief a function that allows you to obtain a value from a sparse matrix
 *         with a binary hint matrix at column col and row row.
 *  This function allows you to retrieve a value from the sparse matrix 
 *  that has a binary hint (index, eye) if a value at a certain postion is 0 or
 *  not. This speeds up retrieval if 0 are called for too.  
 *  \param s_mat sparse matrix to retrieve the value from. 
 *  \param col column of the value to retrieve.
 *  \param row row of the value to retrieve.
 */
double get_value_from_sorted_sparse_matrix_with_eye(sparse_matrix s_mat,
						    int col, int row) {

  int i;
  
  if ( get_value_in_binary_index_at_index(s_mat.eye,
					  (long long)row*(long long)s_mat.range
					  +(long long)col)
       == 1 ) {
    
    i = get_sparse_matrix_index_ordered(s_mat,col,row);
    return(s_mat.values[i]);

  } else {
    return(0.);
  }
  
}

