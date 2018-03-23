/*! \file mesh_from_file.c
 *  \brief This file contains all the code to generate a mesh structure
 *         from an input file. 
 *  The file is used to implement the mesh_from_file() function. Herein all 
 *  the code is found to successfully gerate a structure structre containing
 *  all the necessary to perform the finite element calculation
 */
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hex_fem_solver.h"

/*! \brief helper function to implement getline on systems that do not have it.
 *  Implements GNU getline for systems where it is not available.
 *  The return value is the number of characters on the line read
 *  \param line points to the array containing the line
 *  \param len returns the length of string stored in the array pointed to by 
 *             line
 *  \param fp pointer to the file where to read a line from.
 */
int rc_getline(char **line, size_t *len, FILE *fp)
{
  char *p;
  size_t last = 0;
  
  while (!feof(fp)) {
    if (*line == NULL || last != 0) {
      *len += BUFSIZ;
      *line = realloc(*line, *len);
    }
    p = *line + last;
    memset(p, 0, BUFSIZ);
    if (fgets(p, BUFSIZ, fp) == NULL)
      return -1;
      break;
    last += strlen(p);
    if (last && (*line)[last - 1] == '\n') {
      (*line)[last - 1] = '\0';
      break;
    }
  }
  return last;
}

/*! \brief structure to store a mechanical property before assigning it to 
 *         an element.
 */
typedef struct {
  int id;   /*!< id of the mechanical property*/
  double E;  /*!< Youngs Modulus*/
  double G;  /*!< Shear Modulus*/
  double mu; /*!< Poisson number*/
} mechanical_property;

/*! \brief structure to store a force before assigning it to a node */
typedef struct {
  int id; /*!< node id that this force shall be applied to*/
  int direction; /*!< direction: 0 = x; 1 = y; 2 = z; */
  double f; /*!< value of the force in the given direction*/
} force;

/*! \brief helper function to convert a string to uppercase letters */
char* to_uppercase(char* thestring) {
  int i;
  for(i=0;i<strlen(thestring);i++) {
    thestring[i]=toupper(thestring[i]);
  }
  return thestring;
}

/*! \brief wrapper function arount strtok to avoid segfaults with empty lines *
 *  The my_token function implements a wrapper around strtok that does not 
 *  segfault on empty lines and that converts acquired strings automatically 
 *  to uppercase. Otherwise behavoir is similar to strtok from the c library.
 */
char* my_token(char* init_string, const char* delim, char* buffer) {
  memset(buffer,0,20);
  char* token = strtok(init_string,delim);
  if ( token != NULL ) {
    sscanf(token,"%s",buffer);
    buffer=to_uppercase(buffer);
  } else {
    strncpy(buffer,"ZEROLINE",10);
  }
  return(buffer);
}

/*! \brief function that adds constraints parsed from the input file (i.e.
 *         node fixations to the nodes in question.
 *  The add_constraints_to_nodes() function adds parsed fixations to the
 *  corresponding nodes.
 *  \param nodes pointer to all nodes parsed.
 *  \param n_nodes total number of nodes parsed.
 *  \param fixations array that contains the ids of the nodes to be constrained.
 *  \param n_fixations number of total fixations.
 */
void add_constraints_to_nodes(node* nodes, int n_nodes,
			      int* fixations, int n_fixations) {
  int i,j;
  int node_found;
  for (i=0;i<n_fixations;i++) {
    j=0;
    node_found = 0;
    while (node_found == 0 && n_nodes > j) {
      if(nodes[j].id == fixations[i]) {
	node_found = 1;
	nodes[j].fixed = 1;
      }
      j++;
    }
    if (node_found != 1) {
      printf("Error in input file, constraint %d could not ",i);
      printf("be assigned to \n a node \n");
      exit(-1);
    }
  }
}

/*! \brief function that adds forces parsed from the input file (i.e.
 *         forces applied to nodes to the nodes in question.
 *  The add_forces_to_nodes() function adds parsed fixations to the
 *  corresponding nodes. Note that forces exist in three different directions.
 *  \param nodes pointer to all nodes parsed.
 *  \param n_nodes total number of nodes parsed.
 *  \param forces array that contains the forces and ids of the nodes that
                  these forces shall be applied to.
 *  \param n_forces number of forces.
 */
void add_forces_to_nodes(node* nodes, int n_nodes,
			 force* forces, int n_forces) {
  int i,j,k;
  int node_found;
  for(i=0;i<n_forces;i++) {
    node_found = 0;
    j=0;
    while (node_found == 0 && n_nodes > j ) {
      if (nodes[j].id == forces[i].id) {
	node_found = 1;
	for(k=0;k<3;k++) {
	  if(forces[i].direction == k) {
	    nodes[j].f[k] = forces[i].f;
	  }
	}
      }
      j++;
    }
    if (node_found != 1) {
      printf("Error in input file, force %d could not be assigned to \n",i);
      printf("a node \n");
      exit(-1);
    }
  }
}

/*! \brief helper function that returns the array index for a node with a 
 *         specific id.
 *  The function returns the array index of a node stored in the nodes 
 *  array (first argument) that has the id node_id (second argument). 
 *  The function returns -1 if the node can not be found in the nodes array.
 *  \param nodes array of all the nodes parsed
 *  \param node_id id as specified in the input file
 *  \param n_nodes number of nodes in the nodes array
 */
int node_id_to_index(node* nodes, int node_id, int n_nodes) {
  int i = 0;
  int node_found = 0;

  if ( n_nodes > 0 ) {
    if ( nodes[i].id == node_id ) {
      node_found = 1;
      return 0;
    }
  }

  while (node_found == 0 && n_nodes > i ) {
    if(nodes[i].id == node_id) {
      node_found = 1;
      return i;
    }
    i++;
  }
  // node not found
  return -1;
}

/*! \brief returns the array index of a mechanical property in the 
 *         mechanical_properties array.
 *  This function returns the array index of a mechanical property stored
 *  in the mechanical_properties array (first argument) which has the id
 *  (second argument) as specified in the input file.
 *  \param mechanical_properties an array containing mechanical properties.
 *  \param id the id of the machanical property to be found.
 *  \param number_of_materials the number of mechanical properties stored
 *                             in the mechanical properties array
 */
int find_material_by_id(mechanical_property* mechanical_properties,
			int id,
			size_t number_of_materials) {
  int i = 0;
  int material_found = 0;
  
  if( number_of_materials > 0) {
    if(mechanical_properties[i].id == id) {
      material_found=1;
      return 0;
    }
  }
  
  while(material_found == 0 && number_of_materials > i) {
    i++;
    if(mechanical_properties[i].id == id) {
      material_found=1;
      return i;
    }
  }
  // material not found
  return -1;
}

/*! \brief A helper function to set all displacements and forces in a node
 *         array to zero.
 *  \param nodes a node array.
 *  \param n_nodes the number of nodes in the array.
 */
void initialize_nodes (node* nodes, int n_nodes) {
  int i;

  for(i=0;i<n_nodes;i++) {
    memset(nodes[i].u,0,3); // set displacements to zero
    memset(nodes[i].f,0,3); // set forces to zero
  }
}

/*! \brief The main function that converts a structure defined in an input file 
 *         into into a structure structure.
 *  The read_mesh_from_file() function, builds the structure structure, which
 *  shall hold all values necessary to proceed to a finite element calculation
 *  from the input file.
 *  \param file The input file.
 */
structure read_mesh_from_file(FILE* file) {

  int number_of_nodes = 0;
  int number_of_elements =0;
  int number_of_materials =0;
  int number_of_fixations =0;
  int number_of_forces =0;
  int i,j,k,l,m, ii ;

  int current_material_id;
  
  int nodes_per_element = 8;

  char* line = NULL;
  char buffer[30];
  char* token;
  
  size_t zero_t = 0;
  ssize_t linelength;

  mechanical_property null_material;
  mechanical_property* mechanical_properties;
  int * fixations;
  force* forces;
  structure nodes_and_elements;

  postprocessed_outputs outputs;

  int check;
  
  outputs.out_volumes_before = 0;
  outputs.out_volumes_after = 0;
  outputs.out_displacements = 0;
  outputs.out_forces = 0;
  outputs.out_strain_tensor = 0;
  outputs.out_principal_strains = 0;
  
  null_material.id = 0;
  null_material.E =0;
  null_material.G =0;
  null_material.mu =0;
  
  // initialize counters

  rewind(file);
  
  while ((linelength = rc_getline(&line, &zero_t, file)) != -1) {
    
    // first tokens and outputs
    token = my_token(line, ",",buffer);
    if ( 0 == strncmp(token,"N", 4) ) { number_of_nodes++; }
    if ( 0 == strncmp(token,"E", 4) ) { number_of_elements++; }
    if ( 0 == strncmp(token,"MAT", 4) ) { number_of_materials++; }
    if ( 0 == strncmp(token,"D", 4) ) { number_of_fixations++; }
    if ( 0 == strncmp(token,"F", 4) ) { number_of_forces++; }

    if ( 0 == strncmp(token,"ZOU", 4) ) {
      check = 0;
      token = my_token(NULL,",",buffer);
      if ( 0 == strncmp(token,"VOA",4) ) {
	outputs.out_volumes_after = 1;
	check = 1;
      }
      if ( 0 == strncmp(token,"VOB",4) ) {
	outputs.out_volumes_before = 1;
	check = 1;
      }
      if ( 0 == strncmp(token,"DIS", 4) ) {
	outputs.out_displacements = 1;
	check = 1;
      }
      if ( 0 == strncmp(token,"FOR", 4) ) {
	outputs.out_forces = 1;
	check = 1;
      }
      if ( 0 == strncmp(token,"STE", 4) ) {
	outputs.out_strain_tensor = 1;
	check = 1;
      }
      if ( 0 == strncmp(token,"PST", 4) ) {
	outputs.out_principal_strains = 1;
	check = 1;
      }
      if(!check) {
	printf("Malformed input file ZOU with false specifier \n");
	exit(1);
      }
    }
  }
  nodes_and_elements.nodes =
    (node*)malloc(number_of_nodes*sizeof(node));
  nodes_and_elements.elements =
    (element*)malloc(number_of_elements*sizeof(element));
  mechanical_properties =
    (mechanical_property*)malloc(number_of_materials
				 *sizeof(mechanical_property));
  fixations = (int*) malloc (number_of_fixations*sizeof(int));
  forces = (force*) malloc (number_of_forces*sizeof(force));
  
  if ( nodes_and_elements.nodes == NULL ||
       nodes_and_elements.elements == NULL ) {
    printf("Error allocating memory for nodes and elements \n");
    printf("- read_mesh_from_file() \n");
    exit(1);  
  }

  // parse materials
  
  rewind(file);
  i=0;
  
  while ((linelength = rc_getline(&line, &zero_t, file)) != -1) {
    token = my_token(line,",",buffer);
    if (0 == strncmp(token, "MAT",4)) {
      token = my_token(NULL,",",buffer);
      sscanf(token, "%d", &(mechanical_properties[i].id));
      i++;
    }
  }

  rewind(file);
  
  while ((linelength = rc_getline(&line, &zero_t, file)) != -1) {
    token = my_token(line,",",buffer);
    if (0 == strncmp(token, "MP",4)) {
      token = my_token(NULL,",",buffer);
      if (0 == strncmp(token, "EX",4)) {
	token = my_token(NULL,",",buffer);
	sscanf(token,"%d",&current_material_id);
	i = find_material_by_id(mechanical_properties,
				current_material_id,number_of_materials);
	token = my_token(NULL,",",buffer);
	sscanf(token,"%lf",&(mechanical_properties[i].E));
      }
      if (0 == strncmp(token, "NUXY",4)) {
	token = my_token(NULL,",",buffer);
	sscanf(token,"%d",&current_material_id);
	i = find_material_by_id(mechanical_properties,
				current_material_id,number_of_materials);
	token = my_token(NULL,",",buffer);
	sscanf(token,"%lf",&(mechanical_properties[i].mu));
      }
    }
  }

  // parse nodes and elements
  
  rewind(file);

  // In the loop and only in the loop below 
  // i Nodeindex
  // j Elementindex
  // k constraintindex
  // l current_material_index
  // m current_force_index
  i = j = k = l = m = 0;
  
  while ((linelength = rc_getline(&line, &zero_t, file)) != -1) {

    token = my_token(line,",",buffer);

    // update to current material
    
    if ( 0 == strncmp(token,"MAT",4) ) {
      token = my_token(NULL,",",buffer);
      sscanf(token,"%d",&current_material_id);
      l = find_material_by_id(mechanical_properties,
			      current_material_id, number_of_materials);
    }

    // parse nodes

    if ( 0 == strncmp(token,"N", 4) ) { 
      // element id
      token = my_token(NULL,",",buffer);
      sscanf(token,"%i",&(nodes_and_elements.nodes[i].id));
      // element x,y,z coordinates
      token = my_token(NULL,",",buffer);
      sscanf(token,"%lf",nodes_and_elements.nodes[i].x);
      token = my_token(NULL,",",buffer);
      sscanf(token,"%lf",nodes_and_elements.nodes[i].x+1);
      token = my_token(NULL,"\n",buffer);
      sscanf(token,"%lf",nodes_and_elements.nodes[i].x+2);
      nodes_and_elements.nodes[i].fixed = 0;

      nodes_and_elements.nodes[i].f[0] = 0.;
      nodes_and_elements.nodes[i].f[1] = 0.;
      nodes_and_elements.nodes[i].f[2] = 0.;
      
      i++;
    }

    // parse constraints

    if ( 0 == strncmp(token,"D", 4)) {
      token = my_token(NULL,",",buffer);
      if (1 == sscanf(token,"%d", fixations+k)) {
	k++;
      } else {
	printf("Malformed Input file, constraints unreadable \n");
	exit(-1);
      }
      token = my_token(NULL,",",buffer);
      if ( 0 != strncmp(token,"ALL",4)) {
	printf("Malformed Input file, constraint type not supported \n");
	exit(-1);
      }
      token = my_token(NULL,",",buffer);
      if ( 0 != strncmp(token,"0",4) ) {
	printf("Malformed Input file, constraint type not supported \n");
	exit(-1);
      }
    }
    
    // parse forces

    if ( 0 == strncmp(token,"F", 4)) {
      token = my_token(NULL,",",buffer);
      if (1 == sscanf(token,"%d", &(forces[m].id))) {
	token = my_token(NULL,",",buffer);

	if ( 0 == strncmp(token,"FX",4)) {
	  forces[m].direction = 0;
	}

	if ( 0 == strncmp(token,"FY",4)) {
	  forces[m].direction = 1;
	}

	if ( 0 == strncmp(token,"FZ",4)) {
	  forces[m].direction = 2;
	}
	token = my_token(NULL,",",buffer);
	sscanf(token,"%lf", &(forces[m].f));
	m++;
      } else {
	printf("Malformed Input file, forces unreadable \n");
      }
    }
	
   // parse elements and link material
    
    if ( 0 == strncmp(token,"E", 4) ) {
      nodes_and_elements.elements[j].node_ids =
	(int*)malloc(nodes_per_element*sizeof(int)); 
      
      if ( NULL == nodes_and_elements.elements[j].node_ids ) {
	printf("Could not allocate memory to store nodes for element %d \n",
	       j);
	printf("- read_mesh_from_file() \n");
	exit(1);
      }
      
      for(ii=0; ii< nodes_per_element; ii++) {
	token = my_token(NULL,",",buffer);
	sscanf(token,"%d",nodes_and_elements.elements[j].node_ids+ii);
      }

      nodes_and_elements.elements[j].E_modulus =
	mechanical_properties[l].E;
      nodes_and_elements.elements[j].G_modulus =
	mechanical_properties[l].G;
      nodes_and_elements.elements[j].poisson_number =
	mechanical_properties[l].mu;
      
      j++;
    }
  }
  free(line);
  free(mechanical_properties);

  initialize_nodes(nodes_and_elements.nodes, number_of_nodes);
  
  add_constraints_to_nodes(nodes_and_elements.nodes, number_of_nodes,
			   fixations, number_of_fixations);

  free(fixations);
  
  add_forces_to_nodes(nodes_and_elements.nodes, number_of_nodes,
		      forces, number_of_forces);

  free(forces);
  
  for(i=0;i<number_of_elements;i++) {
    nodes_and_elements.elements[i].node_indicies =
      (int*)malloc(sizeof(int)*8);
    if(nodes_and_elements.elements[i].node_indicies == NULL) {
      printf("Could not allocate Memory for node indicies \n");
      printf("-read_mesh_from_file() \n");
      exit(1);
    }
    for(k=0;k<8;k++) {
      nodes_and_elements.elements[i].node_indicies[k] =
	node_id_to_index(nodes_and_elements.nodes,
			 nodes_and_elements.elements[i].node_ids[k],
			 number_of_nodes);
      if(nodes_and_elements.elements[i].node_indicies[k] == -1) {
	printf("Fatal: Node ID could not be attributed to index \n");
	printf("-read_mesh_form_file() \n");
      }
    }
  }
  nodes_and_elements.n_nodes = number_of_nodes;
  nodes_and_elements.n_elements = number_of_elements;
  nodes_and_elements.outputs = outputs;
  
  printf("nodes %i \n",number_of_nodes);
  printf("elements %i \n",number_of_elements);
  
  return(nodes_and_elements);
}
  
