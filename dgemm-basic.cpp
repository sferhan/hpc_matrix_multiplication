#include "helpers.cpp"

const char* dgemm_desc = "Basic implementation, three-loop dgemm.";

/*
 * This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in column-major format.
 * On exit, A and B maintain their input values.
 */
void square_dgemm(int n, double* A, double* B, double* C) 
{
   // insert your code here: implementation of basic matrix multiplematrix
   for(int i=0; i<n; i++) {
      for(int j=0; j<n; j++) {
         double sum_1r_nc = 0;
         for(int k=0; k<n; k++) {
            sum_1r_nc += (A[col_major_transform(i, k, n)] * B[col_major_transform(k, j, n)]);
         }
         C[col_major_transform(i, j, n)] += sum_1r_nc;
      }
   }
}
