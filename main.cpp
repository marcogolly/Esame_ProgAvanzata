#include <iostream>
#include <chrono>
#include <vector>

#include "parallel_matrix.hpp"
#include "matrix.hpp"

/**
 * @brief Calcola il tempo medio di esecuzione di tutti i metodi della classe Matrix
 * 
 * @param A è la matrice di partenza usata per il calcolo di tutti i metodi
 * @param B è la matrice usata per la molteplicazione in mult_matrix() 
 * @param b è il vettore usato per la molteplicazione in mult_vect() e come vettore dei termini noti in system()
 * @param repetitions è il numero di ripetizioni di ciascuna funzione
 */
template<typename T>
void get_avg_execution_time_Matrix(Matrix<T>& A, const Matrix<T>& B, std::vector<T>& b, const unsigned int repetitions=4) 
{
    auto t0 = std::chrono::high_resolution_clock::now();
    int rank = A.getRank();
    auto t1 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0)/repetitions;
    std::cout << "Tempo di esecuzione di getRank(): " << duration.count() << " microseconds" << std::endl;

    t0 = std::chrono::high_resolution_clock::now();
    int det = A.getDeterminant();
    t1 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0)/repetitions;
    std::cout << "Tempo di esecuzione di getDeterminant(): " << duration.count() << " microseconds" << std::endl;

    t0 = std::chrono::high_resolution_clock::now();
    Matrix<T> prod = A.mult_matrix(B);
    t1 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0)/repetitions;
    std::cout << "Tempo di esecuzione di mult_matrix(): " << duration.count() << " microseconds" << std::endl;

    t0 = std::chrono::high_resolution_clock::now();
    std::vector<T> prod_vect = A.mult_vect(b);
    t1 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0)/repetitions;
    std::cout << "Tempo di esecuzione di mult_vect(): " << duration.count() << " microseconds" << std::endl;

    t0 = std::chrono::high_resolution_clock::now();
    std::vector<T> soluz = A.system(b);
    t1 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0)/repetitions;
    std::cout << "Tempo di esecuzione di system(): " << duration.count() << " microseconds" << std::endl;

    t0 = std::chrono::high_resolution_clock::now();
    Matrix<T> inverse = A.inverse();
    t1 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0)/repetitions;
    std::cout << "Tempo di esecuzione di inverse(): " << duration.count() << " microseconds" << std::endl;
}


/**
 * @brief Calcola il tempo medio di esecuzione di tutte le funzioni della classe Parallel Matrix
 * 
 * @tparam MULTIPLIER è la classe moltiplicatore di matrici considerata dal test
 * @tparam TIMESCALE è l'unità di misura dei tempi calcolati
 * @param A_rows è il numero di righe della prima matrice 
 * @param A_cols è il numero di colonne della prima matrice e il numero di righe della seconda matrice
 * @param B_cols è il numero di colonne della seconda matrice 
 * @param repetitions è il numero di ripetizioni di ciascuna funzione
 */
template<typename MULTIPLIER, typename TIMESCALE=std::chrono::microseconds>
void get_avg_execution_time_ParallelMatrix(const unsigned int A_rows=10, const unsigned int A_cols=10, const unsigned int B_cols=10, const unsigned int repetitions=4)
{
    std::vector<float> A(A_rows*A_cols), B(A_cols*B_cols), b(A_rows);
    std::vector<float> C(A_rows*B_cols), res(A_rows), x(A_cols), inv(A_rows*A_rows);

    srand(time(0)); 
    for (int i = 0; i < A_rows*A_cols; i++) 
        {A[i] = rand()%20;}
    for (int i = 0; i < A_cols*B_cols; i++) 
        {B[i] = rand()%20;}
    for (int i = 0; i < A_rows; i++) 
        {b[i] = rand()%20;}

    int rank = 0;
    float det = 0.0f;
    
    MULTIPLIER multiplier;
    auto t0 = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i<repetitions; ++i) {
        multiplier.getRank(A.data(), A_rows, A_cols, &rank);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Tempo di esecuzione di getRank(): " << std::chrono::duration_cast<TIMESCALE>(t1-t0).count()/repetitions << " microseconds" << std::endl;

    t0 = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i<repetitions; ++i) {
        multiplier.getDeterminant(A.data(), A_rows, A_cols, &det);
    }
    t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Tempo di esecuzione di getDeterminant(): " << std::chrono::duration_cast<TIMESCALE>(t1-t0).count()/repetitions << " microseconds" << std::endl;

    t0 = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i<repetitions; ++i) {
        multiplier.mult_matrix(A.data(), B.data(), C.data(), A_rows, A_cols, B_cols);
    }
    t1 = std::chrono::high_resolution_clock::now();
    std::cout <<  "Tempo di esecuzione di mult_matrix(): " << std::chrono::duration_cast<TIMESCALE>(t1-t0).count()/repetitions << " microseconds" << std::endl;

    t0 = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i<repetitions; ++i) {
        multiplier.mult_vect(A.data(), b.data(), res.data(), A_rows, A_cols);
    }
    t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Tempo di esecuzione di mult_vect(): " << std::chrono::duration_cast<TIMESCALE>(t1-t0).count()/repetitions << " microseconds" << std::endl;

    t0 = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i<repetitions; ++i) {
        multiplier.system(A.data(), b.data(), x.data(), A_rows, A_cols);
    }
    t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Tempo di esecuzione di system(): " << std::chrono::duration_cast<TIMESCALE>(t1-t0).count()/repetitions << " microseconds" << std::endl;

    t0 = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i<repetitions; ++i) {
        multiplier.inverse(A.data(), inv.data(), A_rows, A_cols);
    }
    t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Tempo di esecuzione di inverse(): " << std::chrono::duration_cast<TIMESCALE>(t1-t0).count()/repetitions << " microseconds" << std::endl;
       
}

int main()
{
    
    srand(time(0)); 
    int N = 10;
    Matrix<float> A(N, N),  B(N, N);
    std::vector<float> b(N);
    for (int i=0; i<N; ++i) {
        b[i] = rand()%20;
        for (int j=0; j<N; ++j) {
            A[i][j] = rand()%20;
            B[i][j] = rand()%20;
        }
    }

    std::cout << "Tempi di esecuzione medi della classe " <<  Matrix<float>::name() << ": " << std::endl << std::endl;
    get_avg_execution_time_Matrix(A,B,b);

    std::cout << std::endl<< "---------------------------------------------------------------" << std::endl << std::endl;

    std::cout << "Tempi di esecuzione medi della classe " <<  ParallelMatrix::name() << ": " << std::endl << std::endl;
    get_avg_execution_time_ParallelMatrix<ParallelMatrix>();

    
    return 0;
}
