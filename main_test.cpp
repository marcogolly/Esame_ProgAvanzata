#include <chrono>
#include <functional>
#include <iostream>
#include <vector>
#include "parallel_matrix.hpp"
#include "matrix.hpp"

/**
 * @brief Calcola il tempo medio di esecuzione di tutti i metodi della classe
 * Matrix
 *
 * @param A è la matrice di partenza usata per il calcolo di tutti i metodi
 * @param B è la matrice usata per la molteplicazione in mult_matrix()
 * @param b è il vettore usato per la molteplicazione in mult_vect() e come
 * vettore dei termini noti in system()
 * @param repetitions è il numero di ripetizioni di ciascuna funzione
 */
template <typename T>
void get_avg_execution_time_Matrix(Matrix<T> &A, const Matrix<T> &B,
                                   std::vector<T> &b,
                                   const unsigned int repetitions = 4) {
  std::vector<std::string> nomi_metodi = {"rango", "determinante", "A*B", "A*b", "sistema lineare", "inversa"};

  std::vector<std::function<void()>> method_calls = {
      [&]() {A.getRank(); },      [&]() {A.getDeterminant(); },
      [&]() {A.mult_matrix(B); }, [&]() {A.mult_vect(b); },
      [&]() {A.system(b); },      [&]() {A.inverse(); }};
  
  for (unsigned int i = 0; i < method_calls.size(); i++) {
    auto t0 = std::chrono::high_resolution_clock::now();
    for (unsigned int j = 0; j < repetitions; j++) {
      method_calls[i]();
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    auto duration =std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0) /repetitions;
    std::cout << "Tempo di esecuzione di " << nomi_metodi[i] << ": "
              << duration.count() << " microseconds" << std::endl;
  }
}

template<typename MULTIPLIER, typename TIMESCALE=std::chrono::microseconds>
void get_avg_execution_time_ParallelMatrix(const unsigned int A_rows=10, const unsigned int A_cols=10, const unsigned int B_cols=10, const unsigned int repetitions=4) {
  
  std::vector<float> A(A_rows*A_cols), B(A_cols*B_cols), b(A_rows);
    std::vector<float> C(A_rows*B_cols), res(A_rows), x(A_cols), inv(A_rows*A_rows);

    srand(time(0)); 
    for (int i = 0; i < A_rows*A_cols; i++) 
        {A[i] = rand()%20;}
    for (int i = 0; i < A_cols*B_cols; i++) 
        {B[i] = rand()%20;}
    for (int i = 0; i < A_rows; i++) 
        {b[i] = rand()%20;}

  std::vector<std::string> nomi_metodi = {"rango", "determinante", "A*B", "A*b", "sistema lineare", "inversa"};
  
  MULTIPLIER multiplier;
  std::vector<std::function<void()>> method_calls = {
      [&]() {multiplier.getRank(A.data(), A_rows, A_cols, nullptr); },      [&]() {multiplier.getDeterminant(A.data(), A_rows, A_cols, nullptr); },
      [&]() {multiplier.mult_matrix(A.data(), B.data(), C.data(), A_rows, A_cols, B_cols); }, [&]() {multiplier.mult_vect(A.data(), b.data(), res.data(), A_rows, A_cols);},
      [&]() { multiplier.system(A.data(), b.data(), x.data(), A_rows, A_cols); },      [&]() {multiplier.inverse(A.data(), inv.data(), A_rows, A_cols);}};
  
  for (unsigned int i = 0; i < method_calls.size(); i++) {
    auto t0 = std::chrono::high_resolution_clock::now();
    for (unsigned int j = 0; j < repetitions; j++) {
      method_calls[i]();
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    auto duration =std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0) /repetitions;
    std::cout << "Tempo di esecuzione di " << nomi_metodi[i] << ": "
              << duration.count() << " microseconds" << std::endl;
  }
}

int main() {
  srand(time(0));
  int N = 10;
  Matrix<float> A(N, N), B(N, N);
  std::vector<float> b(N);
  for (int i = 0; i < N; ++i) {
    b[i] = rand() % 20;
    for (int j = 0; j < N; ++j) {
      A[i][j] = rand() % 20;
      B[i][j] = rand() % 20;
    }
  }
  std::cout << "Tempi di esecuzione medi della classe " << Matrix<float>::name()
            << ": " << std::endl
            << std::endl;
  get_avg_execution_time_Matrix(A, B, b);

  std::cout << std::endl
            << "---------------------------------------------------------------"
            << std::endl
            << std::endl;

  //attenzione: per opencl le matrici A,B e il vettore b devono essere float
   std::cout << "Tempi di esecuzione medi della classe " <<
   ParallelMatrix::name() << ": " << std::endl << std::endl;
   get_avg_execution_time_ParallelMatrix<ParallelMatrix>(A, B, b);

  return 0;
}
