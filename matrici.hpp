#include <iostream>
#include <vector>
#include <exception>
#include <algorithm>

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
    int getRows() const 
    { return rows; }

    // funzione che restituisce il numero di colonne (la variabile è privata quindi non è visualizzabile altrimenti)
    int getCols() const 
    { return cols; }

    // funzione che assegna un valore all'elemento alla riga i e colonna j
    void setElement(int i, int j, int val) 
    { matrix[i][j] = val; }

    // funzione che visualizza il valore presente in posizione i,j
    int getElement(int i, int j) const 
    { return matrix[i][j]; }

    //l'istanza dell'oggetto su cui lavoro può essere mofificata
    inline std::vector<T>& operator[](const size_t i)
    {return matrix[i];}

    inline const std::vector<T>& operator[](const size_t i) const
    {return matrix[i];}

    /*
    Matrix gauss(){
        Matrix res(rows,cols);
        res.matrix = matrix;

        for (int i=0; i<rows-1; i++){
            
            if (res[i][i] ==0){
            //flag per vedere se c'è una riga non nulla da scambiare
            bool flag = false;

            for (int l=i;l<rows; l++ )
                if(res[l][i] !=0){
                flag = true;
                //TODO: swap rows i and l
                continue;
                }
            }
            //se non ho trovato una riga da scambiare significa che tutta la colonna i è vuota dalla riga i in giù quindi vado avanti
            if (!flag)
            continue;
            
            for (int j=i+1; j<rows; j++){
            //calcolo lambda dell'algoritmo di gauss per le righe i e j
            float lambda = res[j][i] / res[i][i]; 

            for(int k =i;k<cols; k++){
                res[j][k] -= lambda * res[i][k];
            }
            }
        }
        return res;
    }
    */

    //serve verificare che la matrice sia quadrata altrimenti non si può calcolare il determinante
    


    // funzione per visualizzare la matrice
    void print() const {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                std::cout << matrix[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }


    /*Il rango di una matrice è il numero massimo di elementi indipendenti (cioè non linearmente dipendenti) nella matrice. In altre parole, il rango di una matrice è il numero massimo di colonne o righe indipendenti nella matrice.
Ecco come puoi calcolare il rango di una matrice in C++:
Trasforma la matrice in forma di matrice riflessa (o matrice di scarto). Per fare ciò, puoi utilizzare l'algoritmo di eliminazione di Gauss.
Conta il numero di elementi diversi da zero nella diagonale principale della matrice riflessa. Questo è il rango della matrice originale.
Ecco un esempio di codice che implementa questa procedura:*/
    int rango() 
    {
    int r = 0;
    for (int i = 0; i < rows && r < cols; i++, r++) {
        int pivot = i;
        for (int j = i + 1; j < rows; j++)
            if (abs(matrix[j][r]) > abs(matrix[pivot][r])) 
                {pivot = j;}
        if (pivot != i) 
            {std::swap(matrix[i], matrix[pivot]);}
        if (matrix[i][r] == 0) 
            { r--; continue; }
        for (int j = i + 1; j < rows; j++) {
            double factor = matrix[j][r] / matrix[i][r];
            for (int k = r; k < cols; k++) 
                {matrix[j][k] -= matrix[i][k] * factor;}
        }
    }
    std::cout << "Il rango è: " << r << std::endl;
    return r;
    }
    /*Il codice sopra utilizza l'algoritmo di eliminazione di Gauss per trasformare la matrice in forma di matrice riflessa, quindi conta il numero di elementi diversi da zero sulla diagonale principale per calcolare il rango della matrice.
Spero che questo ti sia stato d'aiuto! Se hai altre domande, non esitare a chiedere.*/


    double determinant() {
    double det = 1;
    for (int i = 0; i < rows; i++) {
        // Trova l'elemento non nullo più a destra sulla riga i-esima
        int pivot = i;
        while (pivot < rows && matrix[i][pivot] == 0) pivot++;
        if (pivot == rows) return 0;  // La matrice è singolare
        if (pivot != i) det *= -1;  // Inverti il segno del determinante se cambi il pivot
        det *= matrix[i][pivot];  // Moltiplica il determinante per il pivot
        for (int j = 0; j < rows; j++) {
        // Sottrai il coefficiente pivot moltiplicato per l'elemento j-esimo della riga i-esima
            matrix[j][i] = (matrix[j][i] * matrix[i][pivot] - matrix[j][pivot] * matrix[i][i]) / matrix[i][pivot];
        }
    } 
    std::cout << "Il determinante è: " << det << std::endl;
    return det;
    }


    friend Matrix<T> operator*(const Matrix<T>& A, const Matrix<T>& B)
    {
        Matrix<T> C(A.getRows(),B.getCols());

        if (A.getRows() == B.getCols()) {
            for (size_t i=0; i<A.getRows(); ++i) {
                for (size_t j=0; j<B.getCols(); ++j) {
                    float partial{0.0};
                    for (size_t k=0; k<A.getCols(); ++k) {
                        partial += A[i][k] * B[k][j];
                    }
                    C[i][j] = partial;
                }
            }
        }
        else {exit;}
        return C;
    }


    Matrix<T> prodottoMatrici(const Matrix<T>& A, const Matrix<T>& B) 
    {
    // Controllo che il numero di colonne della prima matrice sia uguale al numero di righe della seconda matrice
    //if (A.getRows() != B.getCols()) 
    //{throw invalid_argument("Il numero di colonne della prima matrice deve essere uguale al numero di righe della seconda matrice.");}

    // Creo la matrice risultato con un numero di righe pari al numero di righe della prima matrice e un numero di colonne pari al 
    // numero di colonne della seconda matrice
    Matrix<T> C(A.getRows(),B.getCols());

    // Eseguo il prodotto tra le due matrici
    for (size_t i=0; i<A.getRows(); ++i) {
        for (size_t j=0; j<B.getCols(); ++j) {
            for (size_t k=0; k<A.getCols(); ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
    }




    std::vector<T> prod_vect(const Matrix<T>& A, const std::vector<T>& b)
    {
        std::vector<T> c;
        c.size() = A.getRows();

        //if (A.getCols() == b.size()) 
        //{throw invalid_argument("Il numero di colonne della matrice deve essere uguale al numero di elementi del vettore.");}

        for (size_t i=0; i<A.getRows(); ++i) {
            for (size_t j=0; j<A.getCols(); ++j) {
                        c[i] += A[i][j] * b[j];
            }
        }
        return c;
    }



    friend Matrix<T> operator*(const T n, Matrix<T>& A)
    {
        for (size_t i=0; i<A.getRows(); ++i) {
            for (size_t j=0; j<A.getCols(); ++j) {
                A[i][j] = n*A[i][j];                
            }
        }
        return A;
    }

    //Matrix<T>& operator=(const Matrix<T>& orig);



    friend std::ostream& operator<<(std::ostream& os, const Matrix<T>& A)  
    {
    os << "{";
    for (size_t i=0; i<A.getRows(); ++i) {  
        for (size_t j=0; j<A.getCols(); ++j) {
            if (i>0 | j>0) {
                os << ",";
            }  
            os << A[i][j];   
        }
    }
    os << "}"; // non stampo nello standard output ma nello steam os

    return os;    //restituisce riferimento puntatore
    }

    friend std::ostream& operator<<(std::ostream& os, const std::vector<T>& b)  
    {
    os << "[";
    for (size_t i=0; i<b.size(); ++i) {  
            if (i>0) {
                os << ",";
            }  
            os << b[i];   
    }
    os << "]"; // non stampo nello standard output ma nello steam os

    return os;    //restituisce riferimento puntatore
    }

};


/*
// operatore stampa   (tipo di cout è ostream, classe che rappresenta un flusso di dati in output)
template<typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& A)
// secondo parametro oggetto che voglio stampare in output, const perchè voglio solo stampare l'array, non voglio modificarlo
{
os << "{";
for (size_t i=0; i<A.getRows(); ++i) {  //non posso uare len, uso a.size()
    for (size_t j=0; j<A.getCols(); ++j) {
        if (i>0 | j>0) {
            os << ",";
        }  
        os << A[i][j];   // non posso usare a.values[i], uso operatore parentesi quadra
    }
}
os << "}" << std::endl; // non stampo nello standard output ma nello steam os

return os;    //restituisce riferimento puntatore
}*/




