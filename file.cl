
__kernel 
void rank(__global float* A, const unsigned int rows, 
          const unsigned int cols, __global int* rank) 
{
    int i = get_global_id(0);
    int j = get_global_id(1);

    int pivot = i;
    for (int k = i+1; k < rows; k++) {
        if (A[k * cols + j] > A[pivot * cols + j]) {
            pivot = k;
        }
    }
    if (A[pivot * cols + j] != 0) {
        if (i != pivot) {
            for (int k = 0; k < cols; k++) {
                float tmp = A[i * cols + k];
                A[i * cols + k] = A[pivot * cols + k];
                A[pivot * cols + k] = tmp;
            }
        }
        float pivot_val = A[i * cols + j];
        for (int k = j; k < cols; k++) {
            A[i * cols + k] /= pivot_val;
        }
        for (int k = i+1; k < rows; k++) {
            float factor = A[k * cols + j];
            for (int l = j; l < cols; l++) {
                A[k * cols + l] -= A[i * cols + l] * factor;
            }
        }
    }
    // Count number of non-zero rows
    if (i == 0 && j == cols - 1) {
        int r = 0;
        for (int k = 0; k < rows; k++) {
            int nonZero = 0;
            for (int l = 0; l < cols; l++) {
                if (A[k * cols + l] != 0) {
                    nonZero = 1;
                    break;
                }
            }
            r += nonZero;
        }
        *rank = r;
    }
}



__kernel 
void determinant(__global float* A, const unsigned int rows, 
                 const unsigned int cols, __global float* determinant) 
{   
    int row = get_global_id(0);
    int col = get_global_id(1);

    int N=rows;
    *determinant = 1.0f;

    for (int i = 0; i < N; i++) {
        int pivot = i;

        for (int j = i + 1; j < N; j++) {
            if (row == j && A[j * N + col] > A[pivot * N + col]) {
                pivot = j;
            }
        }

        if (row == i) {
            if (pivot != i) {
                *determinant *= -1.0f;
                for (int j = 0; j < N; j++) {
                    float temp = A[i * N + j];
                    A[i * N + j] = A[pivot * N + j];
                    A[pivot * N + j] = temp;
                }
            }
            for (int j = i + 1; j < N; j++) {
                float mult = A[j * N + col] / A[i * N + col];
                for (int k = col; k < N; k++) {
                    A[j * N + k] -= mult * A[i * N + k];
                }
            }
        }
    }

    if (row == 0 && col == 0) {
        for (int i = 0; i < N; i++) {
            *determinant *= A[i * N + i];
        }
    }
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
    int i = get_global_id(0);
    int j = get_global_id(1);

    if (i < rows && j < cols) {
        // Controlla che il pivot scelto non sia nullo
        if (A[j * cols + j] != 0) {
            for (int k = i + 1; k < rows; k++) {
                // Calcola il fattore per la riga corrente
                float factor = A[k * cols + j] / A[j * cols + j];
                for (int l = j; l <= cols; l++) {
                    // Sottrai il fattore moltiplicato per la riga pivot dalla riga corrente
                    A[k * cols + l] = A[k * cols + l] - factor * A[j * cols + l];
                }
                b[k] = b[k] - factor * b[j];
            }
        }
    }

    if (j == cols-1){
        for(int k=rows-1; k>=0; k--) {
            float sum = 0;
            for(int l=k+1; l<cols; l++) {
                sum += A[k * rows + l] * x[l];
            }
            x[k] = (b[k] - sum) / A[k * rows + k];
        }
    }

}

              

__kernel 
void inverse(__global float* A, __global float* inv, unsigned int N) 
{
  
    int i, j, k;
    float pivot, factor;
    int row = get_global_id(0);
    int col = get_global_id(1);

    // initialize Ainv to identity matrix
    if (row == col)
        inv[row * N + col] = 1.0f;
    else
        inv[row * N + col] = 0.0f;

    // forward elimination
    for (i = 0; i < N; i++) {
        if (i == row) {
            pivot = A[row * N + i];
            for (j = i + 1; j < N; j++) {
                if (j == col) {
                    factor = A[j * N + i] / pivot;
                    for (k = i; k < N; k++) {
                        A[j * N + k] -= factor * A[i * N + k];
                    }
                    for (k = 0; k < N; k++) {
                        inv[j * N + k] -= factor * inv[i * N + k];
                    }
                }
            }
        }
    }

    // backward substitution
    if (row == col) {
        for (j = 0; j < N; j++) {
            for (k = i + 1; k < N; k++) {
                inv[row * N + j] -= A[row * N + k] * inv[k * N + j];
            }
            inv[row * N + j] /= A[row * N + i];
        }
    }
  
}    

    

