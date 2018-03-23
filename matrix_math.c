/*! \file matrix_math.c
 *  \brief This file contains math functions around matricies.
 *         Mainly three by three
 *         matricies. Only dense and small general matrix
 *         methods are found in here.
 */
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

/*! \brief calculates the determinant of a three by three matrix.
 *  This function calculates the determinant of a three by three matrix using 
 *  the rule of sarrus.
 */
double det_three_by_three(double* m) {
  // 11 12 13    0 1 2 
  // 21 22 23 == 3 4 5 
  // 31 32 33    6 7 8 
  
  // 11*22*33+12*23*31+13*21*32 -
  // 31*22*13-32*23*11-33*21*12
  
  return(m[0]*m[4]*m[8]+m[1]*m[5]*m[6]+m[2]*m[3]*m[7]-
	 m[6]*m[4]*m[2]-m[7]*m[5]*m[0]-m[8]*m[3]*m[1]);
}

/*! \brief calculates the eigenvalues of real symmetrix three by three matrix.
 *  This function calculates the eigenvalues of real three by three matrix.
 */
double* create_eigenvalues_three_by_three_real_sym(double* m) {

  // 11 12 13    0 1 2 
  // 21 22 23 == 3 4 5 
  // 31 32 33    6 7 8  

  double * eigenvalues = (double*)malloc(sizeof(double)*3);

  double B[9];

  double phi;
  double r; 
  int test = 0;
  
  double p_a = m[1]*m[1]+m[2]*m[2]+m[5]*m[5];
  double q = (m[0]+m[4]+m[8])/3.;

  double p_b =
    (m[0]-q)*(m[0]-q)
    +(m[4]-q)*(m[4]-q)
    +(m[8]-q)*(m[8]-q)+2*p_a;
  
  double p = sqrt(p_b/6);
  double one_over_p = 1/p;

  B[0] = one_over_p * (m[0] - q);
  B[1] = one_over_p * m[1];
  B[2] = one_over_p * m[2];
  B[3] = one_over_p * m[3];
  B[4] = one_over_p * (m[4] - q);
  B[5] = one_over_p * m[5];
  B[6] = one_over_p * m[6];
  B[7] = one_over_p * m[7];
  B[8] = one_over_p * (m[8] - q);

  r = 0.5*det_three_by_three(B);
  
  if( r <= -1 ) {
    phi = M_PI / 3.;
    test = 1;
  }
  if( r >= 1 ) {
    phi = 0.;
    test = 1;
  }
  if(test == 0) {
    phi = acos(r)/3.;
  }

  eigenvalues[0] = q+2*p*cos(phi);
  eigenvalues[2] = q+2*p*cos(phi+(2*M_PI/3.));
  eigenvalues[1] = 3* q - eigenvalues[0] - eigenvalues[2];

  return(eigenvalues);
}

/*! \brief Function that calculates the invers of a three by three matrix 
 *         multiplied by the derminant of the input matrix. 
 *  This function calculates the invers of a three by three matrix 
 *  multiplied by its derminant. This function is used by 
 *  create_invers_three_by_three().
 */
double* create_invers_three_by_three_multiplied_by_determinant(double*m) {
  
  // 11 12 13    0 1 2 
  // 21 22 23 == 3 4 5 
  // 31 32 33    6 7 8 

  double* in = (double*)malloc(9*sizeof(double));

  if(__builtin_expect(in == NULL, 0)) {
    printf("Unable to allocate Memory for invers matrix \n");
    printf(" - create_invers_three_by_three() \n");
  }
  
  in[0] = m[4]*m[8]-m[7]*m[5];
  in[1] = m[2]*m[7]-m[8]*m[1];
  in[2] = m[1]*m[5]-m[4]*m[2];

  in[3] = m[5]*m[6]-m[8]*m[3];
  in[4] = m[0]*m[8]-m[6]*m[2];
  in[5] = m[2]*m[3]-m[5]*m[0];

  in[6] = m[3]*m[7]-m[6]*m[4];
  in[7] = m[1]*m[6]-m[7]*m[0];
  in[8] = m[0]*m[4]-m[3]*m[1];

  return(in);
} 

/*! \brief Function that calculates the invers of a three by three matrix. 
 *  This function calculates the invers of a three by three matrix
 *  multiplied by its derminant.
 */
double* create_invers_three_by_three(double*m) {

  // 11 12 13    0 1 2 
  // 21 22 23 == 3 4 5 
  // 31 32 33    6 7 8 

  int i;
  
  double det = det_three_by_three(m);

  double *in = create_invers_three_by_three_multiplied_by_determinant(m);

  for(i=0;i<9;i++) {
    in[i] = in[i]/det;
  }

  return(in);
} 
