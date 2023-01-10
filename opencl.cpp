#include <iostream>
#include <chrono>
#include <vector>

#include "opencl.hpp"

template<typename MULTIPLIER, typename TIMESCALE=std::chrono::milliseconds>
void get_avg_execution_time(const unsigned int A_rows=1<<6, const unsigned int A_cols=1<<18, const unsigned int B_cols=1<<6, const unsigned int repetitions=4)
{
    std::vector<float> A(A_rows*A_cols, 2);
    /*int* rank{0};
    auto t0 = std::chrono::high_resolution_clock::now();
    MULTIPLIER multiplier;
    for (auto i = 0; i<repetitions; ++i) {
        multiplier.rank(A.data(), A_rows, A_cols, rank);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Tempo di rank: " << std::chrono::duration_cast<TIMESCALE>(t1 - t0).count()/repetitions << std::endl;*/

    float* determinant{0};
    *determinant = 1.0f;
    MULTIPLIER multiplier;
    auto t0 = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i<repetitions; ++i) {
        multiplier.determinant(A.data(), A_rows, A_cols, determinant);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Tempo di determinant: " << std::chrono::duration_cast<TIMESCALE>(t1 - t0).count()/repetitions << std::endl;

    /*std::vector<float> B(A_cols*B_cols, 3), C(A_rows*B_cols);
    auto t0 = std::chrono::high_resolution_clock::now();
    MULTIPLIER multiplier;
    for (auto i = 0; i<repetitions; ++i) {
        multiplier.mult_matrix(A.data(), B.data(), C.data(), A_rows, A_cols, B_cols);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout <<  "Tempo di mult_matrix: " << std::chrono::duration_cast<TIMESCALE>(t1 - t0).count()/repetitions << std::endl;*/

    /*std::vector<float> c(A_cols, 4), res(A_rows);
    auto t0 = std::chrono::high_resolution_clock::now();
    MULTIPLIER multiplier;
    for (auto i = 0; i<repetitions; ++i) {
        multiplier.mult_vect(A.data(), c.data(), res.data(), A_rows, A_cols);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Tempo di mult_vect: " << std::chrono::duration_cast<TIMESCALE>(t1 - t0).count()/repetitions << std::endl;*/

    /*std::vector<float> b(A_rows, 4), x(A_cols);
    auto t0 = std::chrono::high_resolution_clock::now();
    MULTIPLIER multiplier;
    for (auto i = 0; i<repetitions; ++i) {
        multiplier.system(A.data(), b.data(), x.data(), A_rows, A_cols);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Tempo di system: " << std::chrono::duration_cast<TIMESCALE>(t1 - t0).count()/repetitions << std::endl;*/

    /*std::vector<float> H(A_rows*A_rows, 3), inv(A_rows*A_rows);
    auto t0 = std::chrono::high_resolution_clock::now();
    MULTIPLIER multiplier;
    for (auto i = 0; i<repetitions; ++i) {
        multiplier.inverse(H.data(), inv.data(), A_rows);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Tempo di inverse: " << std::chrono::duration_cast<TIMESCALE>(t1 - t0).count()/repetitions << std::endl;*/

       
}

int main()
{
    get_avg_execution_time<ParallelMatrix>();

    //ParallelMatrix obj;
    
    return 0;
}