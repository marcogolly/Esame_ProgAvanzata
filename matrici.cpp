#include <vector>
#include "matrici.hpp"


int main()
{
    std::vector<int> b{1,2,3};
    Matrix<int> M(3,4);
    Matrix<int> A(4,3);
    Matrix<int> B(3,3);
    for (size_t i=0; i<B.getRows(); ++i) {
        for (size_t j=0; j<B.getCols(); ++j) {
            B[i][j] = i+j;
        }
    }
    //Matrix<double> A(3,3);
    //Matrix<float> B(3,3);
    M.setElement(1,1,3);
    M.setElement(2,1,6);
    M.setElement(2,2,0);
    M.setElement(0,1,4);
    M.setElement(0,2,5);
    M.setElement(0,0,6);
    A.setElement(1,1,5);
    A.setElement(2,1,3);
    A.setElement(2,2,7);
    A.setElement(0,1,1);
    A.setElement(0,2,2);
    A.setElement(0,0,3);
    //M.getElement(2,0) = 1; 
    //M.print(); 
    //A.print();
    //B.print();
    
    // std::cout << M*B;
    // std::cout << M*A;
    // std::cout << 2*B[1][2];
    // std::cout << 2*B;
    /*
    std::vector<double> c = prod_vect(A,b);
    std::cout << c;
    Matrix<int> C = prodottoMatrici(M*A);
    std::cout << C;*/
    std::cout << b;

    //std::cout << "M[0][1] = " << M[0][1] << std::endl;
    //std::cout << M << std::endl;
    
    //M.rango();
    //M.determinant();

    return 0;
}
