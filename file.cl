__kernel
//moltiplicazione tra vettore e colonna
void molt(__global float* A, __global float* B, __global float* C,
                   int rows_A, int cols_A)
{
    int row = get_global_id(0);
  int col = get_global_id(1);
  float res = 0.0f;
  for (int k = 0; k < cols_A; ++k) {
    res += A[row * cols_A + k] * B[k * cols_A + col];
  }
  C[row * cols_A + col] = res;
}

//moltiplicazione tra matrice e vettore colonna
__kernel void mult(__global float* A, __global float* x, __global float* y,
                   int rows_A, int cols_A) {
  int row = get_global_id(0);
  float res = 0.0f;
  for (int k = 0; k < cols_A; ++k) {
    res += A[row * cols_A + k] * x[k];
  }
  y[row] = res;
}

//metodo di gauss

//determinare il rango di una matrice con Gauss
//SPERO SIA GIUSTO
__kernel void rank(__global float* A, int rows, int cols, __global int* rank) {
  int row = get_global_id(0);

  // esegui l'eliminazione di Gauss sulla riga corrente
  for (int i = 0; i < rows - 1; i++) {
    if (A[row * cols + i] == 0) {
      // cerca una riga con A[l][i] non nullo
      for (int l = i + 1; l < rows; l++) {
        if (A[l * cols + i] != 0) {
          // scambia la riga l con la riga row
          for (int k = 0; k < cols; k++) {
            float temp = A[row * cols + k];
            A[row * cols + k] = A[l * cols + k];
            A[l * cols + k] = temp;
          }
          break;
        }
      }
    }

    // calcola lambda dell'algoritmo di Gauss per le righe row e i
    float lambda = A[row * cols + i] / A[i * cols + i];

    // sottrai lambda * riga i dalla riga row
    for (int k = i; k < cols; k++) {
      A[row * cols + k] -= lambda * A[i * cols + k];
    }
  }

  // conta le righe non nulle
  if (A[row * cols + rows - 1] != 0) {
    atomic_inc(rank);
  }
}
