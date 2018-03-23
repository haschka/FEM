/*! \file hex_fem_solver.h
 *  \brief The header file containing all the function definitions as well 
 *         as various structure and other definitons.
 *  The hex_fem_solver.h file contains all definitions for all the functions
 *  structures etc. that are shared between the individual c files.
 *  Clear discriptions for functions  are given in the individual
 *  implementation files.
 */

// in the case this shall one day be extended to other meshtypes

// general structure definition
// a system is build like this:
//
//               structure
//            elements  nodes
//
// nodes have ids in order to reference them in the elements structure
// further they also have flat indicies that correspond to a flat node
// array.

/*! \brief The element structure definition.
 *  The element structure contains all the necessary information regarding a 
 *  a single hexaheadral element.
 */
typedef struct {
  int* node_ids;         /*!< node ids as in the input file element definition*/
  int* node_indicies;    /*!< node array indicies as used in c */
  int* neighbours;       /*!< array indicies of adjent cells (if available)*/
  double E_modulus;            /*!< the E modulus ( if available )*/ 
  double G_modulus;            /*!< the G modulus ( if available )*/ 
  double poisson_number;       /*!< the poisson_number (if available)*/
} element;                    

/*! \brief The node structure definition
 * The node structure definition contains all values associated to one node
 */
typedef struct {
  int id;      /*!< The node id as stated in the input file */
  int fixed;   /*!< 0 this node has three degrees of freedom 1 has none */
  double x[3]; /*!< The initial positions as read from the input file */
  double u[3]; /*!< The displacements of each node */
  double f[3]; /*!< The forces applied to these nodes (fixed ones) */
} node;

/*! \brief The structure that knows what to output.
 *  The postprocessed_outputs structure holds all the data to know 
 *  what kinds of values shall be calculated during post processing
 *  and what kinds of values shall be output at the end of all the 
 *  calculations. The values are boolean 0 is false, else is true.
 */
typedef struct {
  int out_volumes_before; /*!< volumes before deformation */
  int out_volumes_after; /*!< volumes after deformation (take displacements) */
  int out_displacements; /*!< print out displacements */
  int out_forces; /*!< print out forces */
  int out_strain_tensor; /*!< print out the strain field for each element */
  int out_principal_strains; /*!< print out the principal strains for 
                                 each element */
} postprocessed_outputs;

/*! \brief The structure structure combines nodes and elements into the 
 *         whole finite element structure to be solved
 * The structure structre definition combines nodes and elements and outputs
 * into one 
 * single structure. This provides a simple handle to redirect almoast all
 * known values of a finite element system to a specific function.
 */
typedef struct {
  node* nodes;       /*!< array of nodes */
  element* elements; /*!< array of elements */
  int n_nodes;       /*!< number of nodes */
  int n_elements;    /*!< number of elements */
  postprocessed_outputs outputs; /*!< structure telling what to output*/
} structure;


// min_max used for node shuffling
// r_mark if not not found elsewhere
typedef struct {
  int min;
  int max;
} min_max;

/*! \brief the sparse matrix structure.
 *  This matrix implements the structure that handles access to 
 *  a square sparse matrix values.
 */
typedef struct {
  int range;         /*!< range of the square sparse matrix*/
  int non_zeros;     /*!< number of non zero entries in the square sparse m*/
  int triangle_type; /*!< 0 = both triangles, 1 = upper -, 2 lower trianlge */
  int* row_indicies; /*!< row indicies array*/
  int* col_indicies; /*!< column indicies array*/
  double* values;    /*!< values array*/
  int has_eye;       /*!< does the matrix have a binary index */
  int has_node_eye;  /*!< does the matrix have a binary index per node (3x3) */
  char* eye;         /*!< the pointer to the binary eye*/
  char* node_eye;    /*!< the pointer to the node eye*/
} sparse_matrix;

/*! \brief two_sparse_matricies structure used to return two matrices from 
 *         certain functions.
 *  This is a small helper structure in order to provide a function to 
 *  return two sparse matricies.
 */
typedef struct {
  sparse_matrix a;
  sparse_matrix b;
} two_sparse_matricies;

// Input
// -----

/*! \brief read_mesh_from_file reads a mesh from the input file,
 *         see implemention in mesh_from_file.c.
 */
structure read_mesh_from_file(FILE* path);

// Matrix Math
// -----------

/*! \brief det_three_by_three calclates the determinant of a three
 *         by three matrix, see implementation in matrix_math.c.
 */
double det_three_by_three(double* matrix);

/*! \brief create_eigenvalues_three_by_three_real_sym calculates
 *         the eigenvalues of a real symmetric three by three matrix,
 *         see implementation in matrix_math.c
 */
double* create_eigenvalues_three_by_three_real_sym(double* matrix);

/*! \brief create_invers_three_by_three calculates
 *         the inverse of a three by three matrix,
 *         see implementation in matrix_math.c
 */
double* create_invers_three_by_three(double* matrix);

/*! \brief create_invers_three_by_three_multiplied_by_determinant
 *         calculates the inverse of a three by three matrix
 *         multiplied by the determinant of the input matrix,
 *         see implementation in matrix_math.c
 */
double* create_invers_three_by_three_multiplied_by_determinant(double* matrix);

// Shape functions
// ---------------

/*! \brief function to obtain the value of the n-th hexaheadral shape function
 *         at point xi, eta, zeta, see implementation in shape_func.c
 */
double shape_func_n(int n, double xi,double eta,double zeta);

/*! \brief function to obtain the value of the n-th hexaheadral shape function
 *         at point xi, eta, zeta, derived to either xi, eta, or zeta,
 *         see implementation in shape_func.c
 */
double shape_func_n_d_g(int n, int g, double xi, double eta, double zeta);

/*! \brief function to obtain the value of the n-th hexaheadral shape function
 *         at point xi, eta, zeta, derived to either x, y, or z,
 *         see implementation in shape_func.c
 */
double shape_func_n_d_l(int n, int g,
			double* invers, // requires the invers jacobi matrix
			double xi, double eta, double zeta);


/*! \brief function to obtain the value of the n-th hexaheadral bubble function
 *         at point xi, eta, zeta, derived to either xi, eta, or zeta,
 *         see implementation in shape_func.c
 */
double bubble_func_n_d_g(int n, int g, double xi, double eta, double zeta);

/*! \brief function to obtain the value of the n-th hexaheadral bubble function
 *         at point xi, eta, zeta, derived to either x, y, or z,
 *         see implementation in shape_func.c
 */
double bubble_func_n_d_l(int n, int g,
			 double* invers, // requires the invers jacobi matrix
			 double xi, double eta, double zeta);


// Node Shuffling and index matricies
// -----------------------------------

/*! \brief function that mappes each index in from a reduced result vector
 *         to a global result vector 
 *         (reduced as fixed nodes have no displacements),
 *         see implementation in node_shuffle.c.
 */
int* generate_reduced_index_to_index_vector_translation(structure object);
/*! \brief function that generates the index matricies, i.e. where a local
 *         stiffness matrix has to be transfered into the globale stiffness
 *         matrix, see implementation in node_shuffle.c.
 */
char** generate_index_and_reduced_index_matrix(structure object,
					       int* reduced_translation_table);

/*! \brief function that counts the number of non zero elements
 *         in an index matrix node_shuffle.c.
 */
int number_of_non_zero_elements_in_index_matrix(char* index_matrix,
						int range);

// Sparse Matrix functions
// -----------------------

// binary eye functions

/*! \brief a function that allocates a binary (hint, 
 *         index, eye) matrix of range range, see implementation in 
 *	   sparse_matrix.c.
 */
char* create_binary_index_matrix(long long range);

/*! \brief a function that sets a value in a binary index ( hint, eye ) matrix,
 *         see implementation in sparse_matrix.c.
 */
void set_value_in_binary_index_at_index(char* bin_index_matrix,long long index);

/*! \brief A function that retrieves a value from a binary index ( hint, eye )
 *         matrix at a specified index, see implemenation in sparse_matrix.c.
 */
int get_value_in_binary_index_at_index(char* bin_index_matrix, long long index);


// create a sparse matrix from an index matrix

/*! \brief function that gernerates a sparse matrix from an index matrix
 *         see implementation in sparse_matrix.c
 */
sparse_matrix generate_sparse_matrix_from_index_matrix(char* index_matrix,
						       int range);

/*! \brief a function that frees the memory of a sparse matrix 
 *         see implementation in sparse_matrix.c.
 */
void free_sparse_matrix(sparse_matrix s_mat);

/*! \brief A function that orders sparse matricies by its index.
 *         see implementaion in sparse_matrix.c.
 */
void order_sparse_matrix(sparse_matrix s_mat);

/*! \brief a function that generates a binray matrix that hints weather 
 *         a value is stored at a certain location in a sparse matrix or not,
 *         see implementation in sparse_matrix.c
 */
char* create_sparse_eye(sparse_matrix s_mat);

/*! \brief a function that allows you to obtain a value from a sparse matrix
 *         at column col and row row, see implementation in sparse_matrix.c
 */
double get_value_from_sorted_sparse_matrix(sparse_matrix s_mat,
					   int col, int row);

/*! \brief a function that allows you to obtain a value from a sparse matrix
 *         with a binary hint matrix at column col and row row, 
 *         see implementation in sparse_matrix.c
 */
double get_value_from_sorted_sparse_matrix_with_eye(sparse_matrix s_mat,
						    int col, int row);


/*! \brief a function that allows to increment a value in a sparse matrix, 
 *         see implementation in sparse_matrix.c
 */
void add_to_value_in_sparse_matrix(sparse_matrix s_matx,
				   int col,
				   int row,
				   double value);

// Stiffnessmatrix and Solve
// -------------------------


/*! \brief function that generates gaussian quadrature points for a 
 *         hexaheadral element, see implemenation in 
 *         stiffness_matrix_for_element.c.
 */
double* generate_quadrature_points();

/*! \brief creates the jacobi matrix for transformations 
 *         partial(x,y,z)/partial(xi,eta,zeta) for a hexaheadral element,
 *         see implementation in stiffness_matrix_for_element.c.
 */	 
double* create_jacobi(node* local_nodes, double xi, double eta, double zeta);

/*! \brief create the stress strain relation (D) for a hexheadral element, see
 *         implementation in stiffness_matrix_for_element.c. 
 */
double* create_strain_stress(element ele);

/*! \brief create the strain stress relation for a hexheadral element, see
 *         implementation in stiffness_matrix_for_element.c. 
 */
double* create_stress_strain(element ele);

/*! \brief create the pseudo strain displacement relation for a 
 *         hexheadral element (D), see
 *         implementation in stiffness_matrix_for_element.c. 
 */
double* create_pseudo_strain_displacement(node* localnodes,
					  double* jac_null_inv,
					  double jac_null_jac_quot,
					  double xi, double eta, double zeta);

/*! \brief create the strain displacement relation for a 
 *         hexheadral element (B), see
 *         implementation in stiffness_matrix_for_element.c. 
 */
double* create_strain_displacement(node* localnodes,
				   double* jacobi,
				   double xi, double eta, double zeta);

/*! \brief create the pseudo stress displacement relation for a 
 *         hexheadral element (D_G), see
 *         implementation in stiffness_matrix_for_element.c. 
 */
double* create_D_G(double* stress_strain, double* pseudo_strain_displacement);

/*! \brief create the lower right part of the elementry stiffnes matrix,
 *         see implementation in stiffness_matrix_for_element.c. 
 */
double* create_GT_D_G(double* pseudo_strain_dispacement, double* D_G);

/*! \brief create the off diagonal part of the elementry stiffnes matrix,
 *         see implementation in stiffness_matrix_for_element.c. 
 */
double* create_BT_D_G(double* strain_displacement, double* D_G);

/*! \brief create the upper left part of the elementry stiffnes matrix,
 *         see implementation in stiffness_matrix_for_element.c. 
 */
double* create_BT_D_B(double* strain_displacement, double* stress_strain);

/*! \brief create an entire condensated stiffness matrix for an element,
 *         see implementation in stiffness_matrix_for_element.c. 
 */
double* create_stiffness_matrix_for_element(structure elements_and_nodes,
					    element ele);

/*! \brief a function that returns the reduced size (the number of nodes 
 *         that are not fixed.), see implementation in 
 *         global_stiffness_matrix.c
 */
int reduced_size(node* nodes, int n_nodes);

/*! \brief a function that returns the reduced size (the number of nodes 
 *         that are not fixed.), see implementation in 
 *         global_stiffness_matrix.c
 */
two_sparse_matricies create_g_and_r_stiffness_matrix_sparse(structure object);

/*! \brief a function solves for the nodal displacements. See implementation
 *         in cg_avx2_without_gather, cg_sse and cg_scalar.
 */
double * obtain_reduced_displacements_cg(structure object,
					 sparse_matrix
					 reduced_stiffness_matrix);

/*! \brief a function that attributes the forces that are exercised to fixed,
 *         nodes, see implementation in global_stiffness_matrix.c
 */
void forces_from_reduced(structure object,
			 sparse_matrix* global_stiffness_matrix);

/*! \brief a function that attributes the nodal displacements from reduced,
 *         without taking account fixed nodes to the form to the global form,
 *         see implemenation in global_stiffness_matrix.c
 */
void node_displacements_from_reduced(structure object, double* reduced);

// postprocessing
// --------------

/*! \brief This function calls the postprocessing functions as defined 
 *         by required outputs in the input file.
 *         see implementation in post_processing.c*/
void postprocessing(structure object);

// for windows only

#ifdef windows
void bzero(void* pointer, size_t len) {
  memset(pointer,0,len);
}
#endif
