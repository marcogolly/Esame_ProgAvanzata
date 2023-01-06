//kernel per trovare MCD che verrà richiamato in Gauss
__kernel 
void MCD(const int a, const int b, __global int* result) {
    if (b == 0) {
        *result = a;
    } else {
        int r = a % b;
        MCD(b, r, result);
    }
}

//kernel gauss 
__kernel void gauss(__global Matrix<T>* res) {
    int pivot_col = 0;
    for (int pivot_row = 0; pivot_col < rows - 1; ++pivot_row) {
        // flag serve a verificare la validità del pivot
        //se è true è tutto a posto
        bool flag = true;

        //se il pivot è nullo devo controllare che tutti gli elementi della colonna i sotto di lui siano nulli
        if (res[pivot_row][pivot_col] == 0) {
            flag = false;

            for (int i = pivot_row; i < rows; ++i) {

                //se trovo una riga con un pivot valido, posso scambiare le righe pivot_row e i
                if (res[i][pivot_col] != 0) {
                    flag = true;
                    atomic_swap(res[i], res[pivot_row]);

                    //siccome ora il pivot non è nullo posso uscire dal ciclo
                    break;
                }
            }
        }
        // se res[pivot_row][pivot_col]==0 e non ho trovato una riga da scambiare significa che
        // dalla riga del pivot in giù, tutta la colonna pivot_col è vuota, quindi passo alla colonna successiva
        if (!flag) {
            pivot_row--;
            pivot_col++;
            continue;
        }

        for (int i = pivot_col + 1; i < rows; i++) {
            // nel caso delle matrici intere moltiplico la riga i per il pivot in modo tale da poter
            // sottrarre alla riga i un multiplo della riga pivot_row garantendo di usare solo numeri interi
            //(is_integer(res[0]))
            if (std::is_integral<T>::value) {
                for (int j = 0; j < cols; ++j) {
                    // questo algoritmo va spesso in overflow anche senza usare numeri enormi, 
                    // questo è inevitabile in alcuni casi siccome C++ non gestisce gli
                    // overflow nelle operazioni tra interi, perciò lo facciamo manualmente
                    //T max= get_max_representable_value(res[0])
                    T max = std::numeric_limits<T>::max() / res[pivot_row][pivot_col];
                    if (max < 0)
                        max = -max;
                    // se vado in overflow alzo un'eccezione
                    if (res[i][j] > max) {
                        std::cout << "Questo " << res[i][j] << " è maggiore di " << max << std::endl;
                        throw std::overflow_error("Errore: impossible utilizzare gauss con questa matrice, valori troppo grandi. Prova a convertire in un tipo più grande (ad esempio longint oppure float")


            // calcolo lambda dell'algoritmo di gauss
            float lambda = res[i][pivot_col] / res[pivot_row][pivot_col];

            // procediamo con l'eliminazione
            for (int j = pivot_col; j < cols; ++j) {
                res[i][j] -= lambda * res[pivot_row][j];
            }

            // per le matrici di interi, semplifichiamo la riga i se possibile
            if (std::is_integral<T>::value) {
                // trovo il massimo comune divisore
                T max_div = res[i][0];
                for (int j = 1; j < cols; ++j) {
                    max_div = MCD(max_div, res[i][j]);
                }
                // nel caso di righe nulle max_div è 0, in quel caso lo poniamo uguale a 1
                if (max_div == 0) {
                    max_div = 1;
                }
                // semplifico la riga
                for (int j = 1; j < cols; ++j) {
                    res[i][j] /= max_div;
                }
            }
        }
      }
    }  
  }
}  
} 




//determino il rango di una matrice usando Gauss
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


__kernel 
void determinant(__global float* A, int rows, int cols, __global float* determinant) {
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
    for (int j = i; j < cols; j++) {
      A[row * cols + j]
    }
  }
}   

//kernel moltiplicazione tra vettore e colonna
__kernel
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

//kernel moltiplicazione tra matrice e vettore colonna
__kernel
void mult(__global float* A, __global float* x, __global float* y,
                   int rows_A, int cols_A) {
  int row = get_global_id(0);
  float res = 0.0f;
  for (int k = 0; k < cols_A; ++k) {
    res += A[row * cols_A + k] * x[k];
  }
  y[row] = res;
}


//kernel per il sistema lineare


//kernel della matrice inversa 
__kernel 
void invert_matrix(__global float* A, __global float* A_inv, int N) {
  int i = get_global_id(0);
  int j = get_global_id(1);
  if (i < N && j < N) {
    // Inizializzazione della matrice L e U
    float L[N][N];
    float U[N][N];
    for (int k = 0; k < N; k++) {
      L[i][k] = (i == k) ? 1 : 0;
      U[k][j] = (k <= j) ? A[k * N + j] : 0;
    }
    // Calcolo della decomposizione LU
    for (int k = 0; k < N; k++) {
      float m = U[k][i] / U[k][k];
      L[i][k] = m;
      for (int l = 0; l < N; l++) {
        U[i][l] = U[i][l] - m * U[k][l];
      }
    }
    // Risoluzione del sistema lineare associato alla matrice A
    float y[N];
    for (int k = 0; k < N; k++) {
      float s = 0;
      for (int l = 0; l < k; l++) {
        s += L[k][l] * y[l];
      }
      y[k] = (A[k * N + N] - s) / L[k][k];
    }
    for (int k = N - 1; k >= 0; k--) {
      float s = 0;
      for (int l = k + 1; l < N; l++) {
        s += U[k][l] * A_inv[l * N + j];
      }
      A_inv[k * N + j] = (y[k] - s) / U[k][k];
    }
  }
}

