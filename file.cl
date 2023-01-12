/**
* Calcola il rango di una matrice tramite il metodo di Gauss
*/
__kernel 
void rank(__global float* A, const unsigned int rows, 
          const unsigned int cols, __global int* rank) 
{   
    int row = get_global_id(0);
    int col = get_global_id(1);

    for (int i=0; i<rows; ++i) {
        int pivot = i;
        
        // determino la riga pivot in cui si trova l'elemento più grande della colonna corrente
        for (int j=i+1; j<rows; ++j) {
            if (row == j && A[j * cols + col] > A[pivot * cols + col]) {
                pivot = j;
            }
        }

        // confronto l'elemento corrente con l'elemento pivot e scambio le righe se necessario
        if (row == i) {
            if (pivot != i) {
                for (int j=0; j<cols; ++j) {
                    float temp = A[i * cols + j];
                    A[i * cols + j] = A[pivot * cols + j];
                    A[pivot * cols + j] = temp;
                }
            }

            // calcolo lambda ed eseguo la sottrazione tra le righe 
            for (int j=i+1; j<rows; ++j) {
                float lambda = A[j * cols + col] / A[i * cols + col];
                for (int k=col; k<cols; ++k) {
                    A[j * cols + k] -= lambda * A[i * cols + k];
                }
            }
        }
    }

    barrier(CLK_GLOBAL_MEM_FENCE);

    // conto le righe non nulle
    if (row == 0 && col == cols-1) {
        for (int k=0; k<rows; ++k) {
            if (A[k * cols + k] != 0) {
                (*rank)++;
            }
        }
    }

    barrier(CLK_GLOBAL_MEM_FENCE);
}


/**
* Calcola il determinante di una matrice tramite il metodo di Gauss
*/
__kernel 
void determinant(__global float* A, const unsigned int rows, 
                 const unsigned int cols, __global float* determinant) 
{   
    int row = get_global_id(0);
    int col = get_global_id(1);

    // solo per matrici quadrate
    int N=rows;
    *determinant = 1.0f;

    for (int i=0; i<N; ++i) {
        int pivot = i;

        // determino la riga pivot in cui si trova l'elemento più grande della colonna corrente
        for (int j=i+1; j<N; ++j) {
            if (row == j && A[j * N + col] > A[pivot * N + col]) {
                pivot = j;
            }
        }

        // confronto l'elemento corrente con l'elemento pivot e scambio le righe se necessario
        if (row == i) {
            if (pivot != i) {
                *determinant *= -1.0f;
                for (int j=0; j<N; ++j) {
                    float temp = A[i * N + j];
                    A[i * N + j] = A[pivot * N + j];
                    A[pivot * N + j] = temp;
                }
            }

            // calcolo lambda ed eseguo la sottrazione tra le righe
            for (int j=i+1; j<N; ++j) {
                float lambda = A[j * N + col] / A[i * N + col];
                for (int k=col; k<N; ++k) {
                    A[j * N + k] -= lambda * A[i * N + k];
                }
            }
        }
    }

    barrier(CLK_GLOBAL_MEM_FENCE);

    // calcolo il determinante come prodotto degli elementi diagonali
    if (row == 0 && col == 0) {
        for (int i=0; i<N; ++i) {
            *determinant *= A[i * N + i];
        }
    }

    barrier(CLK_GLOBAL_MEM_FENCE);
}


/**
* Esegue la molteplicazione tra due matrici
*/
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

  barrier(CLK_GLOBAL_MEM_FENCE);
}


/**
* Esegue la molteplicazione tra una matrice e un vettore
*/
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

  barrier(CLK_GLOBAL_MEM_FENCE);
}


/**
* Risolve un sistema lineare del tipo Ax=b tramite il metodo di Gauss
*/
__kernel
void system(__global float* A, __global float* b, __global float* x, unsigned int rows, unsigned int cols)
{
    int row = get_global_id(0);
    int col = get_global_id(1);

    // solo per matrici quadrate
    int N=rows;

    for (int i=0; i<N; ++i) {
        int pivot = i;

        // determino la riga pivot in cui si trova l'elemento più grande della colonna corrente
        for (int j=i+1; j<N; ++j) {
            if (row == j && A[j * N + col] > A[pivot * N + col]) {
                pivot = j;
            }
        }

        // confronto l'elemento corrente con l'elemento pivot e scambio le righe se necessario
        if (row == i) {
            if (pivot != i) {
                for (int j=0; j<N; ++j) {
                    float temp = A[i * N + j];
                    A[i * N + j] = A[pivot * N + j];
                    A[pivot * N + j] = temp;
                }
            }

            // calcolo lambda ed eseguo la sottrazione tra le righe
            for (int j=i+1; j<N; ++j) {
                float lambda = A[j * N + col] / A[i * N + col];
                for (int k=col; k<N; ++k) {
                    A[j * N + k] -= lambda * A[i * N + k];
                }
                b[j] -= lambda * b[i];
            }
        }
    }

    barrier(CLK_GLOBAL_MEM_FENCE);

    // risolvo il sistema lineare all'indietro
    if (col == N-1 && row == N-1){
        for(int k=N-1; k>=0; k--) {
            float sum = 0;
            for(int l=k+1; l<N; l++) {
                sum += A[k * N + l] * x[l];
            }
            x[k] = (b[k] - sum) / A[k * N + k];
        }
    }
    
    barrier(CLK_GLOBAL_MEM_FENCE);
}

              
/**
* Calcola l'inversa di una matrice 
*/
__kernel 
void inverse(__global float* A, __global float* inv, unsigned int rows, unsigned int cols) 
{
    int N=rows;
    int i, j, k;
    int row = get_global_id(0);
    int col = get_global_id(1);
    for (k = 0; k < N; k++) {
        if(row < N && col < N){
            float factor = A[row * N + k] / A[k * N + k];
            A[row * N + col] -= factor * A[k * N + col];
            inv[row * N + k] = -factor;
        }
    }

    // calculate the inverse
    for (i = 0; i < N; i++) {
        inv[i * N + i] = 1.0f;
    }
    
    for (k = 0; k < j; k++) {
        inv[row * N + col] -= A[row * N + k] * inv[k * N + col];
    }
    
    for (i = 0; i < N; i++) {
        for (j = i + 1; j < N; j++) {
            for (k = j; k < N; k++) {
                inv[i * N + j] -= A[i * N + k] * inv[k * N + j];
            }
        }
    }
    for (i = 0; i < N; i++) {
        float factor = A[i * N + i];
        for (j = 0; j < N; j++) {
            inv[i * N + j] /= factor;
        }
    }
}    

    

