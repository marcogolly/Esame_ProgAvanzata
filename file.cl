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