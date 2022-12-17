#include <iostream>
#include <vector>
#include <exception>

//Definizione della classe
class Matrix {
private:
  // Numero di righe e colonne della matrice. Sono variabili private in modo tale da non poter essere modificate esternamente
  int rows, cols;
  // array 2D contenente gli elementi della matrice
  std::vector<std::vector<int>> matrix;

public:
  // con il costruttore vogliamo definire la matrice fissando la dimensione (r×c)
  Matrix(int r, int c) {
    rows = r; 
    cols = c;
    //inizializzo la matrice (vuota) composta da r righe e c colonne. In particolare ogni riga sarà un array di elementi
    matrix = std::vector<std::vector<int>>(r, std::vector<int>(c));
  }

  // funzione che restituisce il numero di righe (la variabile è privata quindi non è visualizzabile altrimenti)
  int getRows() const { return rows; }

  // funzione che restituisce il numero di colonne (la variabile è privata quindi non è visualizzabile altrimenti)
  int getCols() const { return cols; }

  // funzione che assegna un valore all'elemento alla riga i e colonna j
  void setElement(int i, int j, int val) { matrix[i][j] = val; }

  // funzione che visualizza il valore presente in posizione i,j
  int getElement(int i, int j) const { return matrix[i][j]; }

  //serve verificare che la matrice sia quadrata altriemnti non si può calcolare il determinante
  
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
