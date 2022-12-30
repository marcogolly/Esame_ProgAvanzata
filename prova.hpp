#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <ostream>
#include <vector>
#include <numeric>

//Definizione della classe
template<typename T>
class Matrix {
private:
    // Numero di righe e colonne della matrice. Sono variabili private in modo tale da non poter essere modificate esternamente
    int rows, cols;
    // array 2D contenente gli elementi della matrice
    std::vector<std::vector<T>> matrix;

public:
    // con il costruttore vogliamo definire la matrice fissando la dimensione (r×c)
    Matrix(int r, int c) {
        rows = r; 
        cols = c;
        //inizializzo la matrice (vuota) composta da r righe e c colonne. In particolare ogni riga sarà un array di elementi
        matrix = std::vector<std::vector<T>>(r, std::vector<T>(c));
    }

    Matrix(const std::vector<std::vector<T>>& A) 
    {
        rows = A.size(); 
        cols = A[0].size();
        // Copia i dati dal vector alla matrice
        matrix.resize(rows);
        for (int i = 0; i < rows; i++) {
                matrix[i].resize(cols);
            for (int j = 0; j < cols; j++) {
                matrix[i][j] = A[i][j];
            }
        }
    }

    // funzione che restituisce il numero di righe (la variabile è privata quindi non è visualizzabile altrimenti)
    inline int getRows() const 
    {return rows; }

    // funzione che restituisce il numero di colonne (la variabile è privata quindi non è visualizzabile altrimenti)
    inline int getCols() const 
    {return cols; }

    //l'istanza dell'oggetto su cui lavoro può essere mofificata
    inline std::vector<T>& operator[](const size_t i)
    {return matrix[i];}

    //l'istanza dell'oggetto su cui lavoro non può essere mofificata
    inline const std::vector<T>& operator[](const size_t i) const
    {return matrix[i];}

    /*Matrix<T> gauss() {
        if (std::is_integral<T>::value) {
            throw std::runtime_error("Tipo non valido"); }

        Matrix<T> res(matrix);

        for (size_t i=0; i<rows-1; ++i) {

            bool flag = true;
            if (res[i][i] == 0) {
                // flag per vedere se c'è una riga dove res[l][i] non è nullo
                flag = false;
                for (size_t l=i; l<rows; ++l) {
                    if (res[l][i] != 0) {
                        flag = true;

                        // scambio la riga l con la riga i
                        std::swap(res[l], res[i]);
                        // esco dal for e posso andare avanti
                        break;
                    }
                }
            }
            // se res[i][i] ==0 e non ho trovato una riga da scambiare significa che
            // tutta la colonna i è vuota dalla riga i in giù quindi salto l'indice
            if (!flag)
                continue;

            for (size_t j=i+1; j<rows; ++j) {
                // calcolo lambda dell'algoritmo di gauss per le righe i e j
                T lambda = res[j][i] / res[i][i];

                for (size_t k=i; k<cols; ++k) {
                    res[j][k] -= lambda * res[i][k];
                }
            }
        }
        res.print();
        return res;
        
    }*/

    Matrix<T> gauss() 
    {
        Matrix<T> res(matrix);

        for (int i = 0; i < rows - 1; i++) {

            bool flag = true;
            if (res[i][i] == 0) {
                // flag per vedere se c'è una riga dove res[l][i] non è nullo
                flag = false;
                for (int l = i; l < rows; l++) {
                    if (res[l][i] != 0) {
                    flag = true;

                    // scambio la riga l con la riga i
                    std::swap(res[l], res[i]);
                    // esco dal for e posso andare avanti
                    break;
                    }
                }
            }
            // se res[i][i] ==0 e non ho trovato una riga da scambiare significa che
            // tutta la colonna i è vuota dalla riga i in giù quindi salto l'indice
            if (!flag)
            continue;

            for (int j = i + 1; j < rows; j++) {
                
                if (std::is_integral<T>::value) {
                    for (int k = 0; k < cols; k++)  {  
                        res[j][k] *= res[i][i];
                    }
                }
                // calcolo lambda dell'algoritmo di gauss per le righe i e j
                float lambda = res[j][i] / res[i][i];
                for (int k = i; k < cols; k++) {
                    res[j][k] -= lambda * res[i][k]; 
                }   
            }
        }
        std::cout << res << std::endl;
        return res;
        }

    
    // funzione che calcola il rango della matrice
    int getRank()
    {
        // iniziamo usando l'eliminazione di gauss
        Matrix<T> res = gauss();
        std::cout<< "matrice res: " << std::endl;
        res.print();
        // contiamo le righe non nulle
        int rank = 0;
        for (size_t i=0; i<rows; ++i) {
            bool null_row = false;
            for (size_t j=0; j<cols; ++j) {
                if (res[i][j] != 0) {
                null_row = true;
                break;
                }
            }
        if (null_row) {rank++;}
        }
        std::cout << "Il rango è: " << rank << std::endl;
        return rank;

    }

    T getDeterminant() 
    {
        /*T determinant = 1;
        // Check if the matrix is square
        if (rows != cols) {
            //throw std::invalid_argument("matrice non quadrata");
            std::cout<<"Errore: matrice non quadrata" << std::endl;
        }
        else {
            // uso gauss
            Matrix<float> res = gauss();

            // il determinante è uguale al prodotto dei pivot
            for (size_t i=0; i<rows; ++i) {
                T detNA *= res[i][i];
                
            }
            std::cout << "Il determinante è: " << determinant << std::endl;
        }
        return determinant;*/

        T det = 0;

        if (rows != cols) {
            throw std::invalid_argument("La matrice non è quadrata");
        }
        else {
            // caso base: la matrice è 1x1
            if (rows == 1) {
                return matrix[0][0];
            }
            // caso base: la matrice è 2x2
            if (rows == 2) {
                return matrix[0][0]*matrix[1][1] - matrix[1][0]*matrix[0][1];
            }
            // caso ricorsivo: calcola il determinante utilizzando il metodo di Laplace
            else {
                for (int i = 0; i < cols; i++) {
                    // calcolo il sottodeterminante utilizzando il metodo di Laplace sulla sottomatrice eliminando la riga 0 e la colonna i
                    Matrix<T> submatrix(rows-1, cols-1);
                    for (int j = 1; j < rows; j++) {
                        int submatrix_col = 0;
                        for (int k = 0; k < cols; k++) {
                            if (k == i) continue;
                            submatrix[j - 1][submatrix_col] = matrix[j][k];
                            submatrix_col++;
                        }
                    }
                    T subdet = submatrix.getDeterminant();

                    // aggiungo o sottraggo il sottodeterminante al risultato finale in base al segno dell'elemento corrente (i, 0)
                    if (i % 2 == 0) {
                        det += matrix[0][i] * subdet;
                    } else {
                        det -= matrix[0][i] * subdet;
                    }
                }
            }
        }
        return det;
    }
    
    

    // moltiplicazione riga - colonna matrice matrice
    Matrix<T> mult(const Matrix<T> &B) {
        Matrix<T> res = Matrix(rows, B.cols);

        if (cols != B.rows) {
            throw std::invalid_argument("Errore: il numero di colonne della prima matrice deve essere uguale al numero di righe della seconda matrice");
        } 
        else {
            for (size_t i = 0; i < rows; ++i) {
                for (size_t j = 0; j < B.cols; ++j) {
                    for (size_t k = 0; k < rows; ++k) {
                        res[i][j] += matrix[i][k] * B[k][j];
                    }
                }
            }
        }
        return res;
    }

    // moltiplicazione matrice - vettore colonna
    std::vector<T> mult_vect(const std::vector<T> &b) {
        std::vector<T> res(rows);

        if (cols != b.size()) {
            throw std::invalid_argument("Errore: il numero di colonne della matrice deve essere uguale al numero di elementi del vettore");
        } 
        else {
            for (size_t i = 0; i < rows; ++i) {
                for (size_t j = 0; j < cols; ++j) {
                res[i] += matrix[i][j] * b[j];
                }
            }
        }
        return res;
    }

    // Funzione per aggiungere un vettore colonna alla matrice
    void AddColumn(std::vector<T> &column) {
        // Verifica che il vettore abbia lo stesso numero di elementi delle righe
        // della matrice
        if (column.size() != rows) {
            throw std::invalid_argument("Errore: il vettore colonna ha un numero di elementi diverso dal numero di righe della matrice.");
        } 
        else {
            // Incrementa il numero di colonne della matrice
            cols++;

            // Aggiungi il vettore colonna alla matrice
            for (size_t i = 0; i < rows; i++) {
                matrix[i].push_back(column[i]);
            }
        }
    }

    std::vector<T> RemoveLastColumn() {
        std::vector<T> column(rows);
        for (size_t i = 0; i < rows; ++i) {
            column[i] = matrix[i][rows];
            matrix[i].pop_back(); // metodo pop_back per rimuovere l'ultimo elemento di ogni riga della matrice
        }
        cols = cols - 1;
        
        return column;
    }

    // risoluzione sistema lineare Ax=b    (dim x = num col A)

    // funzione per l'eliminazione di Gauss
    std::vector<T> system(std::vector<T> &b) {

        // copio la matrice di partenza in orig
        Matrix<T> orig(matrix);
        // inizializzo il vettore delle soluzioni 
        std::vector<T> x(cols);

        // controllo che il numero di elementi di b sia uguale al numero di righe della matrice
        if (b.size() != rows) {
            throw std::invalid_argument("Errore: il vettore b deve avere lo stesso numero di elementi del numero di righe della matrice");
        }

        // controllo che la soluzione sia unica
        if (orig.getDeterminant() == 0 | orig.rows!=orig.cols) {
            throw std::invalid_argument("Stop: le soluzioni non sono indipendenti");
        }

        if (std::is_integral<T>::value) {
            throw std::invalid_argument("Tipo non valido");
        }

        else {
            orig.AddColumn(b);
            Matrix<T> S = orig.gauss();            
            std::vector<T> c = S.RemoveLastColumn(); // vettore dei termini noti dopo gauss

            // risolvo il sistema all'indietro
            for (int i = rows - 1; i >= 0; --i) {
            T sum = 0;
            for (int j = i + 1; j < cols; j++) {
                sum += S[i][j] * x[j];
            }
            x[i] = (c[i] - sum) / S[i][i];
            }

            // stampo il vettore delle soluzioni
            std::cout << "Il vettore delle soluzioni è: " << x << std::endl;
            
        }
        return x;
    }


    // matrice inversa (controllo che A*A^-1=matrice identità)
    
    Matrix<T> inverse() 
    {
        Matrix<T> inv(rows, cols);
        Matrix<T> orig(matrix);   // copia della matrice

        if (rows != cols | orig.getDeterminant() == 0) {
            throw std::invalid_argument("Errore: la matrice non è invertibile");
        }

        if (std::is_integral<T>::value) {
            throw std::invalid_argument("Tipo non valido");
        }
        
        else {
            // matrirce identità
            Matrix<T> id(rows, cols);
            for (int i=0; i<rows; ++i) {
                for (int j=0; j<cols; ++j) {
                    if (i==j) 
                        {id[i][j] = 1;}
                    else 
                        {id[i][j] = 0;}
                }
            }

            Matrix<T> B(rows, 2*cols);
            for (int i=0; i<rows; i++) {
                for (int j=0; j<cols; j++) {
                    B[i][j] = orig[i][j];
                }
                for (int j=cols; j < 2*cols; j++) {
                    B[i][j] = id[i][j-cols];
                }
            }

            Matrix<T> G = B.gauss();  
            
            for (size_t i=0; i<rows; ++i) {
                for (size_t j=0; j<i; ++j) {
                    double factor = G[j][i] / G[i][i];
                    for (size_t k=0; k<2*cols; ++k) {
                        G[j][k] -= G[i][k] * factor;
                    }
                }
            }

            // Moltiplica ogni riga per il valore inverso del pivot
            for (size_t i=0; i<rows; ++i) {
                double factor = 1.0 / G[i][i];
                for (size_t j = 0; j < 2*rows; j++) {
                    G[i][j] *= factor;
                }
            }

            // Estrai la matrice inversa dalla destra della matrice estesa
            for (size_t i=0; i<rows; ++i) {
                for (size_t j=0; j<cols; ++j) {
                    inv[i][j] = G[i][j+cols];
                }
            }
            inv.print();
        }
        return inv;
    }

    
    friend Matrix<T> operator*(const Matrix<T>& A, const Matrix<T>& B)
    {
        Matrix<T> C(A.rows,B.cols);

        if (A.cols!=B.rows) {
            throw std::invalid_argument("Errore: il numero di colonne della prima matrice deve essere uguale al numero di righe della seconda matrice");
        }
        else {
            for (size_t i=0; i<A.rows; ++i) {
                for (size_t j=0; j<B.cols; ++j) {
                    for (size_t k=0; k<A.rows; ++k) {
                        C[i][j] += A[i][k] * B[k][j];
                    }
                }
            }
            std::cout << "La matrice prodotto è: " << std::endl << C << std::endl;
        }
        return C;
    }


    // funzione per visualizzare la matrice
    void print() const {
        for (size_t i=0; i<rows; ++i) {
            for (size_t j=0; j<cols; ++j) {
                std::cout << matrix[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    friend std::ostream& operator<<(std::ostream& os, const Matrix<T>& A)  
    {
    for (size_t i=0; i<A.rows; ++i) {  
        for (size_t j=0; j<A.cols; ++j) {
            os << A[i][j] << " ";  
        }
        os << std::endl;
    }
    return os;    //restituisce riferimento puntatore
    }

};


template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& b)
{
    os << "[";
    for (size_t i=0; i<b.size(); ++i) {  
        if (i>0) {
        os << ",";
        }  
        os << b[i];   
    }
    os << "]";
    return os;    
}





