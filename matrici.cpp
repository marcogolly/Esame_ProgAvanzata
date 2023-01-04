#include <vector>
#include <numeric>
//#include "matrici.hpp"
#include "prova.hpp"
#include <chrono>

template<typename T>
void measureExecutionTimes(Matrix<T>& M, const Matrix<T>& B, std::vector<T>& b) 
{
    auto t0 = std::chrono::high_resolution_clock::now();
    int rank = M.getRank();
    auto t1 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0);
    std::cout << "Tempo di esecuzione del metodo getRank(): " << duration.count() << " microseconds" << std::endl;

    t0 = std::chrono::high_resolution_clock::now();
    int det = M.getDeterminant();
    t1 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0);
    std::cout << "Tempo di esecuzione del metodo getDeterminant(): " << duration.count() << " microseconds" << std::endl;

    t0 = std::chrono::high_resolution_clock::now();
    Matrix<T> prod = M.mult(B);
    t1 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0);
    std::cout << "Tempo di esecuzione del metodo mult(): " << duration.count() << " microseconds" << std::endl;

    t0 = std::chrono::high_resolution_clock::now();
    std::vector<T> prod_vect = M.mult_vect(b);
    t1 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0);
    std::cout << "Tempo di esecuzione del metodo mult_vect(): " << duration.count() << " microseconds" << std::endl;

    t0 = std::chrono::high_resolution_clock::now();
    std::vector<T> soluz = M.system(b);
    t1 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0);
    std::cout << "Tempo di esecuzione del metodo system(): " << duration.count() << " microseconds" << std::endl;

    t0 = std::chrono::high_resolution_clock::now();
    Matrix<T> inverse = M.inverse();
    t1 = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0);
    std::cout << "Tempo di esecuzione del metodo inverse(): " << duration.count() << " microseconds" << std::endl;
}


int main()
{
    std::vector<float> b={2,3,3};
    std::vector<float> c={1,2,3,4};
    std::vector<std::vector<float>> d={{1,2},{3,4},{5,6}};
    //Matrix<float> U(d);
    //U.print();
    //Matrix<int> AA(3,3);
    //std::cout << AA << std::endl;
    Matrix<int> M(3,4);
    Matrix<float> A(4,3);
    Matrix<int> B(3,3);
    Matrix<float> E(4,4);
    for (int i=0; i<B.getRows(); ++i) {
        for (int j=0; j<B.getCols(); ++j) {
            B[i][j] = i+j;
        }
    }
    Matrix<int> G(4,4);
    for (int i=0; i<G.getRows(); ++i) {
        for (int j=0; j<G.getCols(); ++j) {
            G[i][j] = i+j-2;
        }
    }
    Matrix<int> P(10,10);
    for (int i=0; i<P.getRows(); ++i) {
        if (i<3) {
            for (int j=0; j<P.getCols(); ++j) {
                P[i][j] = -(i+j) + 3;
            }
        } 
        if (4<i<7) {
            for (int j=0; j<6; ++j) {
                P[i][j] = -(i) + 3;
            }
        }
        else {
            for (int j=0; j<9; ++j) {
                P[i][j] = (i-j) + 2;
            }
        }
    }
    /*P.print(); std::cout << P <<  std::endl;
    int pp = P.getDeterminant(); std::cout << pp << std::endl;
    P.getRank();*/
    Matrix<float> V(3,3);
    V[0][0]=5; V[0][1]=2; V[0][2]=4; V[1][0]=-1; V[1][1]=3; V[1][2]=7; V[2][0]=2; V[2][1]=1; V[2][2]=-3;
    /*V.inverse(); std::cout <<  std::endl;
    V.gauss(); std::cout <<  std::endl;
    V.getRank(); std::cout <<  std::endl;
    int det = V.getDeterminant(); std::cout << det << std::endl; std::cout <<  std::endl;*/
    //V.system(b); std::cout <<  std::endl;
    //G.gauss();
    M[1][1]=3; M[2][1]=6; M[2][2]=1; M[0][1]=4; M[2][0]=5; M[0][0]=2; M[1][3]=1;
    A[1][0]=3; A[2][1]=6; A[2][2]=1; A[0][1]=4; A[2][0]=5; A[0][2]=2;
    E[0][0]=1; E[1][0]=2; E[2][0]=1; E[2][1]=1; E[2][2]=4; E[0][1]=5; E[1][2]=-4; E[0][2]=6; E[1][2]=-4; E[1][3]=7; E[0][3]=-1; E[2][3]=3; E[3][1]=1; E[3][0]=2; E[3][2]=4; E[3][3]=6;
/*
        2 4 0 0                 0 1 2
    M=  0 3 0 1             B=  1 2 3
        5 6 1 0                 2 3 4

        0 4 2                   -2 -1 0 1
    A=  3 0 0               G=  -1  0 1 2
        5 6 1                    0  1 2 3
        0 0 0                    1  2 3 4

        1 5 6 -1
        2 0 -4 7
    E=  1 1 4  3
        2 1 4  6
    */
    measureExecutionTimes(E,A,c);
    //Matrix<float> H = E.inverse();
    //std::cout << E*H << std::endl;
    //A.getDeterminant();
    //std::cout << std::endl;
    //std::cout << M*B << std::endl;
    //M.mult(B);
    /*V.getRank();
    B.getRank();
    M.getRank();  // 3
    G.getRank(); // 2
    A.getRank(); // 3
    E.getRank();*/ // 4
    //G.gauss(); std::cout << std::endl;
    //E.gauss(); 
    //E.intGauss();
    //int det = V.getDeterminant(); std::cout << det << std::endl;//0
    //G.getDeterminant();
    //E.getDeterminant();
    
    
    //G.mult(A);  std::cout << std::endl;
    //A.mult(G);  std::cout << std::endl;
    //M.mult(A);  std::cout << std::endl;
    //M.mult(B);  std::cout << std::endl;
    //std::cout << b << std::endl;
    //B.print();
    //A.AddColumn(c);
    //A.print();
    //A.RemoveLastColumn(); A.print();
    //A.print();
    /*A.AddColumn(c);
    A.print();   std::cout << std::endl;
    Matrix<float> S = A.gauss();   std::cout  << std::endl;
    S.print();   std::cout << std::endl;
    std::vector<float> f = S.RemoveLastColumn();
    std::cout << f << std::endl;
    S.print();*/
    //G.system(c);
    //V.system(b); std::cout <<  std::endl;
    //E.system(c);
    //E.inverse();

    //B.system(b);
    //A.print();
    //A.system(c);
    //G.system(c);
    //std::cout << A.system(c);
    //E.inverse();
    //G.inverse();  //condizione det=0
    //A.inverse();
    //B.mult_vect(b); std::cout << std::endl;
    //A.mult_vect(b); std::cout << std::endl;
    //M.mult_vect(c); std::cout << std::endl;
    //G.mult_vect(c); 
    //std::cout << M*c << std::endl;
    //A.print();
    //A.getRank();
    //A.getDeterminant();
    //B.print();
    //B.getRank();
    //B.getDeterminant();
    //M.getElement(2,0) = 1; 
    //M.print(); 
    //A.print();
    //B.print();
    //G.print();
    
    // std::cout << M*B;
    // std::cout << M*A;
    // std::cout << 2*B[1][2];
    // std::cout << 2*B;
    /*
    std::vector<double> c = prod_vect(A,b);
    std::cout << c;
    Matrix<int> C = prodottoMatrici(M*A);
    std::cout << C;*/
    //std::cout << b;

    //std::cout << "M[0][1] = " << M[0][1] << std::endl;
    //std::cout << M << std::endl;
    
    //M.rango();
    //M.determinant();

    /*auto t0 = std::chrono::high_resolution_clock::now();
    M.mult(A);
    auto t1 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t1-t0);
    std::cout << "Execution time: " << duration.count() << " microseconds" << std::endl;*/

    return 0;
}
