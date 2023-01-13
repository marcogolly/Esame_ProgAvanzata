/**
* Calcola il rango di una matrice tramite il metodo di Gauss
*/
__kernel 
void getRank(__global float* A, const unsigned int rows, 
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
void getDeterminant(__global float* A, const unsigned int rows, 
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
void system(__global float* A, __global float* b, __global float* x, 
            unsigned int rows, unsigned int cols)
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
    if (row == N-1 && col == N-1){
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
    int row = get_global_id(0);
    int col = get_global_id(1);

    int N=rows;

    // inizializzo la matrice completa A|I
    __global float* completa = (__global float*) (A);

    // Creiamo una matrice estesa A|I
    for(int i=0; i<N; ++i) {
        // prima metà della matrice completa è uguale a A
        for(int j=0; j<N; ++j) {
            completa[row * 2*N + j] = A[row * N + j]; 
        }
        //seconda metà della matrice completa è la matrice identità
        for(int j=N; j<2*N; ++j) {
            if (j-N == i) {
                completa[row * 2*N + j] = 1; 
            } else {
                completa[row * 2*N + j] = 0;
            }
        }
    }

    for (int i=0; i<N; ++i) {
        int pivot = i;

        // determino la riga pivot in cui si trova l'elemento più grande della colonna corrente
        for (int j=i+1; j<N; ++j) {
            if (row == j && completa[j * N + col] > completa[pivot * N + col]) {
                pivot = j;
            }
        }

        // confronto l'elemento corrente con l'elemento pivot e scambio le righe se necessario
        if (row == i) {
            if (pivot != i) {
                for (int j=0; j<N*2; ++j) {
                    float temp = completa[i * N + j];
                    completa[i * N + j] = completa[pivot * N + j];
                    completa[pivot * N + j] = temp;
                }
            }

            // calcolo lambda ed eseguo la sottrazione tra le righe
            for (int j=i+1; j<N; ++j) {
                float lambda = completa[j * N + col] / completa[i * N + col];
                for (int k=col; k<N*2; ++k) {
                    completa[j * N + k] -= lambda * completa[i * N + k];
                }
            }
        }
    }

    barrier(CLK_GLOBAL_MEM_FENCE);

    // annullo elementi sopra la diagonale in modo da ottenere 
    // la matrice identità nella parte sinistra della matrice completa
    if (row == N-1 && col == N-1) {
        for(int k=N-1; k>=0; k--) {
            if(row == k) {
                float lambda = completa[k * N + k];
                for(int j=k; j<N*2; ++j) {
                    completa[k * N + j] /= lambda;
                }
                for(int i=0; i<k; ++i) {
                    if(row == i) {
                        lambda = completa[i * N + k];
                        for(int j=k; j<N*2; ++j) {
                            completa[i * N + j] -= lambda * completa[k * N + j];
                        }
                    }
                }
            }
        }
    }

    barrier(CLK_GLOBAL_MEM_FENCE);

    // salva solo la parte destra della matrice completa A|I nella matrice inv
    for(int j=N; j<N*2; j++) {
        inv[row * N + j - N] = completa[row * N + j];
    }

    barrier(CLK_GLOBAL_MEM_FENCE);
}


    

