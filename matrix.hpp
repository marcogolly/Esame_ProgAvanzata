#ifndef __matrix_hpp__  
#define __matrix_hpp__ 

#include <iostream>
#include <ostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <exception>


/**
 * @brief Calcola il massimo comun divisore tra due numeri interi con l'algortimo di Euclide
 * 
 * @param a è il primo numero intero per il calcolo del MCD
 * @param b è il secondo numero intero per il calcolo del MCD
 * @return il massimo comun divisore tra a e b
*/
int MCD(const int a, const int b) 
{
    if (b==0)
        return a;
    return MCD(b, a%b);
}

/**
 * @brief Una classe per la rappresentazione delle matrici
 *
 * @param T è il tipo degli elementi della matrice
*/
template<typename T>
class Matrix {
private:
    uint rows; //! è il numero di righe della matrice
    uint cols; //! è il numero di colonne della matrice
    std::vector<std::vector<T>> matrix; //! è l'array 2D contenente gli elementi della matrice

public:
    /**
     * @brief Un costruttore: costruisce matrice nulla di dimensione rxc ( ogni riga sarà un array di elementi)
    */    
    Matrix(int r, int c) {
        rows = r; 
        cols = c;
        if (r == 0 || c == 0) {
            throw std::invalid_argument("Matrice non valida");
        }
        matrix = std::vector<std::vector<T>>(r, std::vector<T>(c));
    }

    /**
     * @brief Un costruttore di copia dato un vector di vector
     * 
     * @param orig è l'istanza originale della copia
    */
    Matrix(const std::vector<std::vector<T>>& orig) 
    {
        rows = orig.size(); 
        cols = orig[0].size();
        if (rows==0 || cols==0) {
            throw std::invalid_argument("Matrice non valida");
        }
        // Copia i dati dal vector alla matrice
        matrix.resize(rows);
        for (int i=0; i<rows; i++) {
                matrix[i].resize(cols);
            for (int j=0; j<cols; j++) {
                matrix[i][j] = orig[i][j];
            }
        }
    }

    /**
     * @brief Un costruttore di copia data una matrice
     * 
     * @param orig è l'istanza originale della copia
    */
    Matrix(const Matrix<T>& orig) 
    {
        rows = orig.rows;
        cols = orig.cols;
        if (rows==0 || cols==0) {
            throw std::invalid_argument("Matrice non valida");
        }
        matrix = orig.matrix;
    }

    /**
     * Ottienne il numero di righe
     * 
     * @return il numero di righe della matrice
    */
    inline const uint getRows() const 
    {return rows; }

    /**
     * Ottienne il numero di colonne
     * 
     * @return il numero di colonne della matrice
    */
    inline const uint getCols() const 
    {return cols; }

    /**
     * @brief Accede all'i-esimo elemento della matrice
     * 
     * @param i è l'indice dell'elemento ambito
     * @return un riferiemento all'elemento della matrice
     */    
    inline std::vector<T>& operator[](const size_t i)
    {return matrix[i];}

    /**
     * @brief Accede all'i-esimo elemento della matrice
     * 
     * @param i è l'indice dell'elemento ambito
     * @return un riferiemento all'elemento della matrice
     */    
    inline const std::vector<T>& operator[](const size_t i) const
    {return matrix[i];}

    /**
     * @brief Calcola la matrice ridotta di una matrice con il metodo dell'eliminazione di Gauss
     * 
     * @return la matrice ridotta
    */
    Matrix<T> gauss(){   
        // usiamo una copia della matrice in modo tale da non modificare l'originale
        Matrix<T> res(*this);
        int pivot_col=0;
        for (int pivot_row=0; pivot_col<rows-1; ++pivot_row) {
            // flag serve a verificare la validità del pivot
            //se è true è tutto a posto
            bool flag = true;

            //se il pivot è nullo devo controllare che tutti gli elementi della colonna i sotto di lui siano nulli
            if (res[pivot_row][pivot_col]==0) {
                flag = false;

                for (int i=pivot_row; i<rows; ++i) {
                
                    //se trovo una riga con un pivot valido, posso scambiare le righe pivot_row e i
                    if (res[i][pivot_col] != 0) {
                        flag = true;
                        std::swap(res[i], res[pivot_row]);

                        //siccome ora il pivot non è nullo posso uscire dal ciclo
                        break;
                    }
                }
            }
            // se res[pivot_row][pivot_col]==0 e non ho trovato una riga da scambiare significa che
            // dalla riga del pivot in giù, tutta la colonna pivot_col è vuota, quindi passo alla colonna successiva
            if (!flag){
                pivot_row--;
                pivot_col++;
                continue;
            }
            
            for (int i = pivot_col + 1; i < rows; i++) {
                // nel caso delle matrici intere moltiplico la riga i per il pivot in modo tale da poter
                // sottrarre alla riga i un multiplo della riga pivot_row garantendo di usare solo numeri interi
                if (std::is_integral<T>::value) {
                    for (int j=0; j<cols; ++j) {
                        // questo algoritmo va spesso in overflow anche senza usare numeri enormi, 
                        // questo è inevitabile in alcuni casi siccome C++ non gestisce gli
                        // overflow nelle operazioni tra interi, perciò lo facciamo manualmente
                        T max = std::numeric_limits<T>::max()/res[pivot_row][pivot_col];
                        if (max<0)
                            max = -max;
                        // se vado in overflow alzo un'eccezione
                        if (res[i][j]>max){
                            std::cout << "Questo " << res[i][j] << " è maggiore di " << max <<std::endl;
                            throw std::overflow_error("Errore: impossible utilizzare gauss con questa matrice, valori troppo grandi. Prova a convertire in un tipo più grande (ad esempio longint oppure float");
                        }
                        res[i][j] *= res[pivot_row][pivot_col];
                    }
                }

                // calcolo lambda dell'algoritmo di gauss
                float lambda = res[i][pivot_col] / res[pivot_row][pivot_col];

                // procediamo con l'eliminazione
                for (int j=pivot_col; j<cols; ++j) {
                    res[i][j] -= lambda * res[pivot_row][j];
                }

                // per le matrici di interi, semplifichiamo la riga i se possibile
                if (std::is_integral<T>::value) {
                    // trovo il massimo comune divisore
                    T max_div = res[i][0];
                    for (int j=1; j<cols; ++j) {
                        max_div = MCD(max_div, res[i][j]);
                    }
                    // nel caso di righe nulle max_div è 0, in quel caso lo poniamo uguale a 1
                    if (max_div==0) {
                        max_div = 1;
                    }
                    // semplifico la riga
                    for (int j=1; j<cols; ++j) {
                        res[i][j] /= max_div;
                    }
                }
            }
            pivot_col++;
        }
        return res;
    }

    
    /**
     * @brief Calcola il rango della matrice tramite l'eliminazione di gauss
     *
     * @return il rango della matrice
    */  
    const int getRank() 
    {
        Matrix<T> res = Matrix(gauss());
        
        // contiamo le righe non nulle
        int rank = 0;
        for (uint i=0; i<rows; ++i) {
            bool null_row = true;
            for (uint j=0; j<cols; ++j) {
                if (res[i][j] != 0) {
                null_row = false;
                break;
                }
            }
        if (!null_row) {rank++;}
        }
        return rank;
    }

    /**
     * @brief Calcola il determinante di una matrice con il metodo di Laplace
     * 
     * @return il determinante della matrice
    */
    const T getDeterminant() 
    {
        T det = 0;

        if (rows != cols) {
            throw std::invalid_argument("La matrice non è quadrata");
        }
        else {
            // caso base: la matrice è 1x1
            if (rows==1) {
                return matrix[0][0];
            }
            // caso base: la matrice è 2x2
            if (rows==2) {
                return matrix[0][0]*matrix[1][1] - matrix[1][0]*matrix[0][1];
            }
            // caso ricorsivo: calcola il determinante utilizzando il metodo di Laplace
            else {
                for (uint i=0; i<cols; ++i) {
                    // calcolo il sottodeterminante utilizzando il metodo di Laplace sulla sottomatrice eliminando la riga 0 e la colonna i
                    Matrix<T> submatrix(rows-1,cols-1);
                    for (uint j=1; j<rows; ++j) {
                        int submatrix_col = 0;
                        for (uint k=0; k<cols; ++k) {
                            if (k==i) 
                                continue;
                            submatrix[j-1][submatrix_col] = matrix[j][k];
                            submatrix_col++;
                        }
                    }
                    T subdet = submatrix.getDeterminant();

                    // aggiungo o sottraggo il sottodeterminante al risultato finale in base al segno dell'elemento corrente (i, 0)
                    if (i%2==0) {
                        det += matrix[0][i] * subdet;
                    } else {
                        det -= matrix[0][i] * subdet;
                    }
                }
            }
        }
        return det;
    }
    

    /**
     * @brief Moltiplicazione tra due matrici righe per colonne
     * 
     * @param B è la matrice usata per la molteplicazione
     * @return la matrice risultato della molteplicazione
    */
    Matrix<T> mult_matrix(const Matrix<T> &B) const
    {
        Matrix<T> res = Matrix(rows, B.cols);

        if (B.cols==0 || B.rows==0) {
            throw std::invalid_argument("Matrice B nulla");
        }
        if (cols != B.rows) {
            throw std::invalid_argument("Errore: il numero di colonne della prima matrice deve essere uguale al numero di righe della seconda matrice");
        } 
        else {
            for (uint i=0; i<rows; ++i) {
                for (uint j=0; j<B.cols; ++j) {
                    for (uint k=0; k<rows; ++k) {
                        res[i][j] += matrix[i][k] * B[k][j];

                    }
                }
            }
        }
        return res;
    }

    /**
     * @brief Moltiplicazione tra una matrice e un vettore
     * 
     * @param b è il vettore usato per la molteplicazione
     * @return il vettore risultato della molteplicazione
    */
    std::vector<T> mult_vect(const std::vector<T> &b) const
    {
        std::vector<T> res(rows);

        if (cols != b.size()) {
            throw std::invalid_argument("Errore: il numero di colonne della matrice deve essere uguale al numero di elementi del vettore");
        } 
        else {
            for (uint i=0; i<rows; ++i) {
                for (uint j=0; j<cols; ++j) {
                res[i] += matrix[i][j] * b[j];
                }
            }
        }
        return res;
    }

    /**
     * @brief Aggiunge un vettore colonna ad una matrice
     * 
     * @param column è il vettore da aggiungere alla matrice
    */
    void AddColumn(std::vector<T> &column) 
    {
        // verifico che il vettore abbia lo stesso numero di elementi delle righe della matrice
        if (column.size() != rows) {
            throw std::invalid_argument("Errore: il vettore colonna ha un numero di elementi diverso dal numero di righe della matrice.");
        } 
        else {
            // incremento il numero di colonne della matrice
            cols++;

            // aggiungo il vettore colonna alla matrice
            for (int i=0; i<rows; ++i) {
                matrix[i].push_back(column[i]);
            }
        }
    }

    /**
     * @brief Rimuove l'ultima colonna di una matrice
     * 
     * @return la colonna rimossa
    */
    std::vector<T> RemoveLastColumn() 
    {
        std::vector<T> column(rows);
        for (int i=0; i<rows; ++i) {
            column[i] = matrix[i][rows];
            matrix[i].pop_back(); // metodo pop_back per rimuovere l'ultimo elemento di ogni riga della matrice
        }
        // aggiorno il numero di colonne
        cols = cols-1;
        
        return column;
    }


    /**
     * @brief Risolve un sistema lineare del tipo Ax=b attraverso l'eliminazione di Gauss
     * 
     * @param b è il vettore dei termini noti
     * @return il vettore delle soluzioni
    */
    std::vector<T> system(std::vector<T> &b) const
    {
        // copio la matrice di partenza in orig
        Matrix<T> orig(*this);
        // inizializzo il vettore delle soluzioni 
        std::vector<T> x(cols);

        // controllo che il numero di elementi di b sia uguale al numero di righe della matrice
        if (b.size() != rows) {
            throw std::invalid_argument("Errore: il vettore b deve avere lo stesso numero di elementi del numero di righe della matrice");
        }

        // controllo che la soluzione sia unica
        if (orig.getDeterminant() == 0 || orig.rows!=orig.cols) {
            throw std::invalid_argument("Stop: le soluzioni non sono indipendenti");
        }

        if (std::is_integral<T>::value) {
            throw std::invalid_argument("Impossible risolvere il sistema lineare di interi");
        }

        else {
            // aggiungo il vettore dei termini noti alla matrice per creare la matrice completa
            orig.AddColumn(b);
            // eliminazioen di gauss sulla matrice completa
            Matrix<T> S = orig.gauss();  
            // estraggo il vettore dei termini noti dopo gauss
            std::vector<T> c = S.RemoveLastColumn(); 
            // risolvo il sistema all'indietro
            for (int i=rows-1; i>=0; --i) {
                T sum = 0;
                for (int j=i+1; j<cols; j++) {
                    sum += S[i][j] * x[j];
                }
                x[i] = (c[i] - sum) / S[i][i];
            }
        }
        return x;
    }

    /**
     * @brief Calcola l'inversa di una matrice con il metodo di Gauss-Jordan
     * 
     * @return la matrice inversa
    */
    Matrix<T> inverse() const
    {
        Matrix<T> inv(rows, cols);
        Matrix<T> orig(*this);   // copia della matrice

        if (rows != cols || orig.getDeterminant() == 0) {
            throw std::invalid_argument("Errore: la matrice non è invertibile");
        }

        if (std::is_integral<T>::value) {
            throw std::invalid_argument("Impossible calcolare l'inversa di una matrice di interi");
        }
        
        else {
            // matrirce identità
            Matrix<T> id(rows, cols);
            for (uint i=0; i<rows; ++i) {
                for (uint j=0; j<cols; ++j) {
                    if (i==j) 
                        {id[i][j] = 1;}
                    else 
                        {id[i][j] = 0;}
                }
            }

            // matrice completa
            Matrix<T> B(rows, 2*cols);
            for (uint i=0; i<rows; ++i) {
                // copio la matrice di partenza
                for (uint j=0; j<cols; ++j) {
                    B[i][j] = orig[i][j];
                }
                // copio la matrice d'identità
                for (uint j=cols; j<2*cols; ++j) {
                    B[i][j] = id[i][j-cols];
                }
            }

            Matrix<T> G = B.gauss();  
            
            // sottrai la riga corrente dalle righe sottostanti per annullare i coefficienti sotto la diagonale
            for (uint i=0; i<rows; ++i) {
                for (uint j=0; j<i; ++j) {
                    double factor = G[j][i] / G[i][i];
                    for (uint k=0; k<2*cols; ++k) {
                        G[j][k] -= G[i][k] * factor;
                    }
                }
            }

            // moltiplica ogni riga per il valore inverso del pivot
            for (uint i=0; i<rows; ++i) {
                double factor = 1.0 / G[i][i];
                for (uint j = 0; j < 2*rows; ++j) {
                    G[i][j] *= factor;
                }
            }

            // estrai la matrice inversa dalla destra della matrice completa
            for (uint i=0; i<rows; ++i) {
                for (uint j=0; j<cols; ++j) {
                    inv[i][j] = G[i][j+cols];
                }
            }
        }
        return inv;
    }

    /**
     * @brief Stampa la matrice
    */
    void print() const 
    {
        for (uint i=0; i<rows; ++i) {
            for (uint j=0; j<cols; ++j) {
                std::cout << matrix[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    /**
     * @brief Nome della classe
     * 
     * @return il nome della classe
     */
    static std::string name()
    {
        return "Matrix";
    }

};

/**
 * @brief Stampa matrice su uno stream
 * 
 * @tparam T è il tipo degli elementi della matrice da stampare
 * @param os è lo stream in cui stampare
 * @param A è la matrice da stampare
 * @return il riferimento allo stream in cui abbiamo stampato
*/
template<typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& A)  
{
    os << "{" << std::endl;
    for (uint i=0; i<A.getRows(); ++i) {
        os << " {";
        for (uint j=0; j<A.getCols(); ++j) {
            if (j > 0) {
                os << ", ";
            }
            os << A[i][j];
        }
        os << "}" << std::endl; 
    }
    os << "}";
    return os;   
}

/**
 * @brief Stampa vettore su uno stream
 * 
 * @tparam T è il tipo degli elementi del vettore da stampare
 * @param os è lo stream in cui stampare
 * @param b è il vettore da stampare
 * @return il riferimento allo stream in cui abbiamo stampato
*/
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& b)
{
    os << "[";
    for (const auto& value : b) {
    os << value;
        if (&value != &b.back()) {
        os << ",";
        }
    }
    os << "]";
    return os;    
}


#endif  // __matrix_hpp__
