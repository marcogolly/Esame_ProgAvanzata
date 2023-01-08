
void AddColumn(__global float* A, __global float* column,  
               unsigned int rows, unsigned int cols) 
{
  int i = get_global_id(0);
  A[i * (cols+1) + cols] = column[i];
  cols += 1;
}

void RemoveLastColumn(__global float* A, __global float* column, unsigned int rows, unsigned int *cols) 
{
  int i = get_global_id(0);  
  column[i] = A[i * (*cols) + (*cols-1)];
  *cols -= 1;
}


__global float* RemoveColumn(__global float* A, unsigned int rows, unsigned int cols) 
{
  __global float* column;
  for (int i = 0; i < rows; i++) {
      column[i] = A[i * (cols) + (cols)];
      A[i * (cols) + (cols)] = 0;
  }
  (cols)--;
  return column;
}

void gauss(__global float* A, unsigned int rows, unsigned int cols) 
{
  int pivot_row = get_global_id(0);
    int pivot_col = get_global_id(1);
    if (pivot_row >= rows - 1 || pivot_col >= rows) return;

  bool flag = true;
  if (A[pivot_row * cols + pivot_col] == 0) {
    flag = false;
    for (int i = pivot_row; i < rows; ++i) {
      if (A[i * cols + pivot_col] != 0) {
        flag = true;
        for (int j = pivot_col; j < cols; ++j) {
          float temp = A[i * cols + j];
          A[i * cols + j] = A[pivot_row * cols + j];
          A[pivot_row * cols + j] = temp;
        }
        break;
      }
    }
  }

  if (!flag) {
    pivot_row--;
    pivot_col++;
    return;
  }

  for (int i = pivot_col + 1; i < rows; ++i) {
    float lambda = A[i * cols + pivot_col] / A[pivot_row * cols + pivot_col];
    for (int j = pivot_col; j < cols; ++j) {
      A[i * cols + j] -= lambda * A[pivot_row * cols + j];
    }
  }
}


__kernel 
void rank(__global const float* A, const unsigned int rows, 
          const unsigned int cols, __global int* rank) 
{
  int i = get_global_id(0);
  int j = get_global_id(1);

  gauss(A,rows,cols);

  if (A[i*cols + j] != 0) {
    (*rank)++;
  }
}



__kernel 
void determinant(__global const float* A, const unsigned int rows, 
                 const unsigned int cols, __global float* determinant) 
{
  int i = get_global_id(0);

  gauss(A,rows,cols);

  *determinant *= A[i*cols + i];
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

    __global float* column = RemoveColumn(A, rows, cols);
    
    for (int i=rows-1; i>=0; --i) {
      float sum = 0;
      for (int j=i+1; j<cols; j++) {
          sum += A[i * cols + j] * x[j];
      }
      x[i] = (column[i] - sum) / A[i * cols + i];
    }      
}



__kernel 
void inverse(__global float* A, __global float* inv, unsigned int N) 
{
  int row = get_global_id(0);
  int col = get_global_id(1);

  int i = get_global_id(0);
  int j = get_global_id(1);

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

  if (row<N && j<N) {
    inv[row*N + col] = B[row*N + col+N];
  }
}    

    

