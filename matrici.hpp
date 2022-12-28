#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <ostream>
#include <vector>

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


    // restituisce la matrice convertita in float
    Matrix<float> toFloat() 
    {
        Matrix<float> res(rows, cols);

        for (size_t i=0; i<rows; ++i) {
            for (size_t j=0; j<cols; ++j) {
                res[i][j] = static_cast<float>(matrix[i][j]);
            }
        }
        return res;
    }


    Matrix<float> gauss() 
    {
        Matrix<float> res(rows, cols);
        res = toFloat();

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
                float lambda = res[j][i] / res[i][i];

                for (size_t k=i; k<cols; ++k) {
                    res[j][k] -= lambda * res[i][k];
                }
            }
        }
        res.print();
        return res;
    }

    
    Matrix<float> Gauss1() {
        // verifica che la matrice sia quadrata
        //if (numeroRighe != numeroColonne) {throw std::invalid_argument("La matrice non è quadrata.");}

        Matrix<float> res(rows,cols);
        res = toFloat();
        // eliminazione di Gauss
        for (size_t i=0 ; i<rows; i++) {
            // cerca il pivot massimo nella colonna i-esima
            float pivot = i;
            for (size_t j=i+1; j<cols; j++) {
                if (std::abs(res[j][i]) > std::abs(res[pivot][i])) {
                    pivot = j;
                }
            }

            if (pivot != i) {
                for (size_t j=0; j<cols; j++) {
                // scambia la riga pivot con la riga i-esima
                std::swap(res[i][j], res[pivot][j]);
                }
            }

            // sottrai la riga i-esima dalle righe j-esime per j > i
            for (size_t j=0; j<rows; j++) {
                auto m = res[j][i] / res[i][i];
                for (size_t k=0; k<cols; k++) {
                    res[j][k] -= m * res[i][k];
                }
            }
        }
        return res;
    }    
    
    

    int getRank()
    {
        // Iniziamo usando l'eliminazione di gauss
        Matrix<float> res = gauss();
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
        T determinant = 1;
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
                determinant *= res[i][i];
            }
            std::cout << "Il determinante è: " << determinant << std::endl;
        }
        return determinant;
    }


    //moltiplicazione riga - colonna matrice matrice
    Matrix<T> mult(const Matrix<T>& B)
    {
        Matrix<T> res = Matrix(rows, B.cols);

        if (cols!=B.rows) {
            std::cout << "Errore: il numero di colonne della prima matrice deve essere uguale al numero di righe della seconda matrice" << std::endl;
        }
        else {
            for (size_t  i=0; i<rows; ++i){
                for (size_t j=0; j<B.cols; ++j){
                    for (size_t k =0; k< rows; ++k){
                        //std::cout<< matrix[i][k] * B[k][j]<< "\n";
                        res[i][j] += matrix[i][k] * B[k][j];
                    }
                }
            }
            std::cout << "La matrice prodotto è: " << std::endl << res << std::endl;
        }
        return res;
    }


    //moltiplicazione matrice - vettore colonna
    std::vector<T> mult_vect(const std::vector<T>& b)
    {
        std::vector<T> res(rows);

        if (cols!=b.size()) {
            std::cout << "Errore: il numero di colonne della matrice deve essere uguale al numero di elementi del vettore" << std::endl;
        }
        else {
            for (size_t  i=0; i<rows; ++i){
                for (size_t j=0; j<cols; ++j){
                    res[i] += matrix[i][j] * b[j];
                }
            }
            std::cout << "Il vettore prodotto è: " << std::endl << res << std::endl;
        }
        return res;
    }


    // Funzione per aggiungere un vettore colonna alla matrice
    void AddColumn(std::vector<float>& column) 
    {
        // Verifica che il vettore abbia lo stesso numero di elementi delle righe della matrice
        if (column.size() != rows) {
            std::cout << "Errore: il vettore colonna ha un numero di elementi diverso dal numero di righe della matrice." << std::endl;
        }
        else {
            // Incrementa il numero di colonne della matrice
            cols++;

            // Aggiungi il vettore colonna alla matrice
            for (size_t i=0; i<rows; i++) {
                matrix[i].push_back(column[i]);
            }
        }
    }

    std::vector<T> RemoveLastColumn() 
    {
        std::vector<T> column(rows);
        for (size_t i=0; i<rows; ++i) {
            column[i] = matrix[i][rows];   
            matrix[i].pop_back();  // metodo pop_back per rimuovere l'ultimo elemento di ogni riga della matrice
            
        }
        rows = rows;
        cols = cols-1;
        //std::cout << column << std::endl;
        //std::cout << matrix;
        return column;
    }

    
    

    // risoluzione sistema lineare Ax=b    (dim x = num col A)

    // funzione per l'eliminazione di Gauss
    std::vector<float> system(std::vector<float>& b)  
    {
        std::vector<float> x(cols);

        // Verifica che il numero di elementi di b sia uguale al numero di righe della matrice
        if (b.size() != rows) {
            //throw std::invalid_argument("Errore: il vettore b deve avere lo stesso numero di elementi del numero di righe della matrice");
            std::cout << "Errore: il vettore b deve avere lo stesso numero di elementi del numero di righe della matrice" <<std::endl;
        }

        else {
            // copio la matrice di partenza in orig
            Matrix<T> orig(rows, cols);
            for (size_t i=0; i<rows; ++i) {
                for (size_t j=0; j<cols; ++j) {
                    orig[i][j] = matrix[i][j];
                }
            }
            orig.toFloat();
            orig.AddColumn(b);
            Matrix<float> res = orig.gauss();
            Matrix<float> S(rows, cols);  // copio res in S
            for (size_t i=0; i<rows; ++i) {
                for (size_t j=0; j<cols; ++j) {
                    S[i][j] = res[i][j];
                }
            }
            
            std::vector<float> c = res.RemoveLastColumn();  // termini noti dopo gauss  
            
            if (res.getDeterminant() == 0) {
                std::cout << "Stop: le soluzioni non sono indipendenti" << std::endl;
            }

            else {
                std::cout << "Avanti: " << std::endl;
                // risolvo il sistema all'indietro
                for (int i=rows-1; i>=0; --i) {
                    double sum = 0;
                    for (int j=i+1; j<cols; j++) {
                        sum += S[i][j] * x[j];
                    }
                    x[i] = (c[i] - sum) / S[i][i];
                }

                // stampo il vettore delle soluzioni
                //std::cout << "Il vettore delle soluzioni è: " << x << std::endl;
                for (int i=0; i<cols; ++i) {
                    std::cout << "Il vettore delle soluzioni è: " << x[i] << std::endl; 
                }
            }
        }
        return x;    
    }


    // calcola la trasposta della matrice (utile per calcolare l'inversa)
    Matrix transpose() 
    {
        Matrix res(cols, rows);
        for (size_t i=0; i<rows; ++i) {
            for (size_t j=0; j<cols; ++j) {
                res[j][i] = matrix[i][j];
            }
        }
        return res;
    }


    // matrice inversa (controllo che A*A^-1=matrice identità)
    
    Matrix<float> inverse() 
    {
        Matrix<float> inv(rows, cols);

        // matrirce identità
        Matrix<T> id(rows, cols);
        for (int i=0; i<rows; ++i) {
            for (int j=0; j<cols; ++j) {
                if (i==j) {id[i][j] = 1;}
                else {id[i][j] = 0;}
            }
        }
        //id.print();  std::cout << std::endl;

        Matrix<T> B(rows, 2*cols);
        for (int i=0; i<rows; i++) {
            for (int j=0; j<cols; j++) {
                B[i][j] = matrix[i][j];
            }
            for (int j=cols; j < 2*cols; j++) {
                B[i][j] = id[i][j-cols];
            }
        }

        if (rows != cols) {// | A.getDeterminant() == 0) {
            //throw std::invalid_argument("matrice non invertibile");
            std::cout<<"Errore: la matrice non è invertibile" << std::endl; }
        
        else {
            B.print(); std::cout << std::endl;
            Matrix<float> G = B.gauss();  std::cout << std::endl;
            
            for (size_t i=0; i<rows; ++i) {
                // Scegli il pivot
                /*int pivot = i;
                for (int j = i+1; j < rows; j++) {
                    if (fabs(G[j][i]) > fabs(G[pivot][i])) {
                        pivot = j;
                    }
                }
                if (pivot != i) {
                // Scambia le righe pivot e i
                    for (int j = 0; j < 2*cols; j++) {
                        std::swap(G[i][j], G[pivot][j]);
                    }
                }*/
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

    
    /*friend Matrix<T> operator=(const Matrix<T>& A)
    {
        for (int i=0; i<rows; i++) {
            for (int j=0; j<cols; j++) {
                A[i][j] = matrix[i][j];
            }
    }*/
    
    
    friend Matrix<T> operator*(const Matrix<T>& A, const Matrix<T>& B)
    {
        Matrix<T> C(A.rows,B.cols);

        if (A.cols!=B.rows) {
            std::cout << "Errore: il numero di colonne della prima matrice deve essere uguale al numero di righe della seconda matrice" << std::endl;
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



    /*
    friend Matrix<T> operator*(const T n, Matrix<T>& A)
    {
        for (size_t i=0; i<A.getRows(); ++i) {
            for (size_t j=0; j<A.getCols(); ++j) {
                A[i][j] = n*A[i][j];                
            }
        }
        return A;
    }
    */


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

std::ostream& operator<<(std::ostream& os, const std::vector<float>& b)
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



