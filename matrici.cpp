#include <vector>
#include "matrici.hpp"


int main()
{
    std::vector<float> b={2,3,3};
    std::vector<float> c={1,2,3,4};
    Matrix<int> M(3,4);
    Matrix<int> A(4,3);
    Matrix<int> B(3,3);
    Matrix<float> E(4,4);
    for (size_t i=0; i<B.getRows(); ++i) {
        for (size_t j=0; j<B.getCols(); ++j) {
            B[i][j] = i+j;
        }
    }
    Matrix<int> G(4,4);
    for (size_t i=0; i<G.getRows(); ++i) {
        for (size_t j=0; j<G.getCols(); ++j) {
            G[i][j] = i+j-2;
        }
    }
    //Matrix<double> A(3,3);
    //Matrix<float> B(3,3);
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
    //Matrix<float> H = E.inverse();
    //std::cout << E*H << std::endl;
    //M.getRank();  // 3
    //B.getDeterminant();
    //std::cout << std::endl;
    //std::cout << M*B << std::endl;
    //M.mult(B);
    //G.getRank(); // 2
    //A.getRank();
    //G.getDeterminant();  //0
    //A.getDeterminant();
    //E.getDeterminant();
    //E.getRank();
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
    //E.system(c);
    //M.system(b);
    //A.print();
    //A.system(c);
    G.system(c);
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

    return 0;
}
