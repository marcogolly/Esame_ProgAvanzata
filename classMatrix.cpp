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
  int getRows() const { return rows; }

  // funzione che restituisce il numero di colonne (la variabile è privata quindi non è visualizzabile altrimenti)
  int getCols() const { return cols; }

  
  
  //override dell'operatore [][] che permette di scrivere Matrix[i][j]
  inline std::vector<T>& operator[](size_t i) {
    return matrix[i];
  }

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
  
  
  // funzione per visualizzare la matrice
  void print() const {
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        std::cout << matrix[i][j] << " ";
      }
      std::cout << std::endl;
    }
  }
};



int main(){
  srand(time(0));  // Seed per il random

    Matrix<float> m(3,3);
    
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            //se voglio inizializzare una matrice random
            //m[i][j] = rand();

            //inizializzazione di prova
            m[i][j] = i+1+j;
        }
    }
    std::cout<<"matrice m: \n";
    m.print();
    std::cout<<"matrice gauss: \n";
    m.gauss().print();

    m.print();

}
