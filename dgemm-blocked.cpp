#include "helpers.cpp"

const char* dgemm_desc = "Blocked dgemm.";

void load_block_in_cache(int n, int block_size, int br, int bc, double* X) {
   double* temp_mem = new double(block_size);
   for(int i=0; i<block_size; i++) {
      memcpy(temp_mem, X+((br*block_size+i)*sizeof(double)), block_size*sizeof(double));
   }
}

void block_operate(
   int n, int xr, int xc, int yr, int yc, int blocks_count,
   int block_size, double* X, double* Y, double* Z
) {
   for(int i=0; i<block_size; i++) {
      for(int j=0; j<block_size; j++) {
         double sum_1r_nc = 0;
         for(int k=0; k<block_size; k++) {
            sum_1r_nc += (X[col_major_transform(xr*block_size+i, xc*block_size+k, n)] * Y[col_major_transform(yr*block_size+k, yc*block_size+j, n)]);
         }
         Z[col_major_transform(xr*block_size+i, yc*block_size+j, n)] += sum_1r_nc;
      }
   }
}

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm_blocked(int n, int block_size, double* A, double* B, double* C) 
{
   // insert your code here
   int nb = n/block_size;

   for(int i=0; i<nb; i++) {
      for(int j=0; j<nb; j++) {
         load_block_in_cache(n, block_size, i, j, C);
         for(int k=0; k<nb; k++) {
            load_block_in_cache(n, block_size, i, k, A);
            load_block_in_cache(n, block_size, k, j, B);
            block_operate(n, i, k, k, j, nb, block_size, A, B, C);
         }
      }
   }
}
