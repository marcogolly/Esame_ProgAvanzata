void AddColumn(__global float* A, __global float* column,  
               unsigned int rows, unsigned int cols) 
{
  int i = get_global_id(0);
  A[i * (cols+1) + cols] = column[i];
  cols += 1;
}


void gauss(__global float* A, unsigned int rows, unsigned int cols) 
{
  int pivot_row = get_global_id(0);
  int pivot_col = get_global_id(1);

  if (pivot_row >= rows || pivot_col >= cols) {
    return;
  }

  // Scegli il pivot di questa riga
  for (int i = pivot_col + 1; i < cols; i++) {
    if (fabs(A[pivot_row * cols + i]) > fabs(A[pivot_row * cols + pivot_col])) {
      pivot_col = i;
    }
  }

  // Scambia questa riga con la riga che contiene il pivot più grande
  if (pivot_col != pivot_row) {
    for (int i = 0; i < cols; i++) {
      float temp = A[pivot_row * cols + i];
      A[pivot_row * cols + i] = A[pivot_col * cols + i];
      A[pivot_col * cols + i] = temp;
    }
  }

  // Rendi zero gli elementi sottostanti il pivot
  for (int i = pivot_row + 1; i < rows; i++) {
    float lambda = A[i * cols + pivot_col] / A[pivot_row * cols + pivot_col];
    for (int j = 0; j < cols; j++) {
      A[i * cols + j] -= lambda * A[pivot_row * cols + j];
    }
  }
}


__kernel 
void rank(__global float* A, const unsigned int rows, 
          const unsigned int cols, __global int* rank) 
{
  int i = get_global_id(0);
  int j = get_global_id(1);

  gauss(A,rows,cols);

  if (A[i*cols + j] != 0) {
    *rank += 1;
  }
}



__kernel 
void determinant(__global float* A, const unsigned int rows, 
                 const unsigned int cols, __global float* determinant) 
{
    int i = get_global_id(0);
    int j = get_global_id(1);
    int pivot_col = 0;
    int pivot_row = 0;
    float result = 1.0;

    for (pivot_row = 0; pivot_col < cols-1; ++pivot_row) {
        // trova il massimo pivot
        int max_row = pivot_row;
        for (i = pivot_row+1; i < rows; ++i) {
            if (fabs(A[i * rows + pivot_col]) > fabs(A[max_row * rows + pivot_col])) {
                max_row = i;
            }
        }
        //se trovo una riga con un pivot valido, posso scambiare le righe pivot_row e max_row
        if (max_row != pivot_row) {
            for (j = pivot_col; j < cols; ++j) {
                float temp = A[pivot_row * rows + j];
                A[pivot_row * rows + j] = A[max_row * rows + j];
                A[max_row * rows + j] = temp;
            }
            result = -result;
        }
        // procediamo con l'eliminazione
        float pivot = A[pivot_row * rows + pivot_col];
        for (i = pivot_row + 1; i < rows; i++) {
            float lambda = A[i * rows + pivot_col] / pivot;
            for (j = pivot_col; j < cols; j++) {
                A[i * rows + j] = A[i * rows + j] - lambda * A[pivot_row * rows + j];
            }
        }
        pivot_col++;
    }

    //calcoliamo il determinante come prodotto dei diagonal principal
    for(i = 0; i < rows; i++)
    {
        result *= A[i * rows + i];
    }

    *determinant = result;

}   




__kernel 
void mult_matrix(__global const float *A, __global const float *B, __global float *C,
                 const unsigned int A_rows, const unsigned int A_cols, const unsigned int B_cols)
{
  int i = get_global_id(0);
  int j = get_global_id(1);

  float partial = 0;
  for (uint k=0; k<A_cols; ++k) {
    partial += A[i*A_cols + k] * B[k*B_cols + j];
  }

  C[i*B_cols + j] = partial;
}



__kernel 
void mult_vect(__global const float *A, __global const float *b, __global float *res,
               const unsigned int rows, const unsigned int cols)
{
  int i = get_global_id(0);

  float partial = 0;
  for (uint j=0; j<cols; ++j) {
    partial += A[i*cols + j] * b[j];
  }

  res[i] = partial;
}


__kernel
void system(__global float* A, __global float* b, __global float* x, unsigned int rows, unsigned int cols)
{
  int row = get_global_id(0);
  int col = get_global_id(1); 

  AddColumn(A, b, rows, cols);

  gauss(A, rows, cols);  
  
  float sum = 0;
  for (int j = row + 1; j < cols - 1; j++) {
    sum += A[row * cols + j] * x[j];
  }
  x[row] = ((row*cols + (cols-1)) - sum) / A[row*cols + col];    
}

              

__kernel 
void inverse(__global float* A, __global float* inv, unsigned int N) 
{
  int row = get_global_id(0);
  int col = get_global_id(1);

  __global float* id;
  if (row == col) {
    id[row*N + col] = 1.0f;
  } else {
    id[row*N + col] = 0.0f;
  }

  __global float* B;
  B[row*(2*N) + col] = A[row*N + col];
  B[row*(2*N) + col + N] = id[row*N + col];

  gauss(B, N, 2*N);

  if (row<N && col<row) {
    double factor = B[col*N + row] / B[row*N + row];
    for (uint k=0; k<2*N; ++k) {
      B[col*N + k] -= B[row*N + k] * factor;
    }
  }

  if (row<N) {
    double factor = 1.0 / B[row*N + row];
    for (uint j=0; j<2*N; ++j) {
      B[row*N + j] *= factor;
    }
  }

  if (row<N && col<N) {
    inv[row*N + col] = B[row*N + col+N];
  }
}    

    

