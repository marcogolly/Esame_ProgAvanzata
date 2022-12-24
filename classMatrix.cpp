#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <ostream>
#include <vector>

// Definizione della classe
template <typename T> class Matrix {
private:
  // Numero di righe e colonne della matrice. Sono variabili private in modo
  // tale da non poter essere modificate esternamente
  int rows, cols;
  // array 2D contenente gli elementi della matrice
  std::vector<std::vector<T>> matrix;

public:
  // con il costruttore vogliamo definire la matrice fissando la dimensione
  // (r×c)
  Matrix(int r, int c) {
    rows = r;
    cols = c;
    // inizializzo la matrice (vuota) composta da r righe e c colonne. In
    // particolare ogni riga sarà un array di elementi
    matrix = std::vector<std::vector<T>>(r, std::vector<T>(c));
  }

  // restituisce la matrice convertita in float
  Matrix<float> toFloat() {
    Matrix<float> res(rows, cols);

    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        res[i][j] = static_cast<float>(matrix[i][j]);
      }
    }
    return res;
  }

  // funzione che restituisce il numero di righe (la variabile è privata quindi
  // non è visualizzabile altrimenti)
  int getRows() const { return rows; }

  // funzione che restituisce il numero di colonne (la variabile è privata
  // quindi non è visualizzabile altrimenti)
  int getCols() const { return cols; }

  // override dell'operatore [][] che permette di scrivere Matrix[i][j]
  inline std::vector<T> &operator[](size_t i) { return matrix[i]; }

  // funzione che restituisce la matrice dopo l'eliminazione di gauss
  // nota: la funzione restituisce una matrice float anche per le matrici di
  // interi
  Matrix<T> intGauss() {
    Matrix<T> res(rows, cols);
    res =*this;

    for (int i = 0; i < rows - 1; i++) {

      bool flag = true;
      std::cout<<"res ii "<<res[i][i]<<" "<<i<<std::endl<<res;
      if (res[i][i] == 0) {
        // flag per vedere se c'è una riga dove res[l][i] non è nullo
        flag = false;
        for (int l = i; l < rows; l++) {
          if (res[l][i] != 0) {
            flag = true;

            // scambio la riga l con la riga i
            std::cout<<res;
            std::swap(res[l], res[i]);
            std::cout<<res;
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
        if (res[j][i] % res[i][i] != 0)
          for (int k = 0; k < cols; k++)
            res[j][k] *= res[i][i];
        // calcolo lambda dell'algoritmo di gauss per le righe i e j
        T lambda = res[j][i] / res[i][i];

        for (int k = i; k < cols; k++) {
          res[j][k] -= lambda * res[i][k];
        }
      }
    }
    return res;
  }

  Matrix<T> gauss() {
    if (std::is_same<T, int>::value) {
      return intGauss();
    }

    Matrix<T> res(rows, cols);
    res = *this;

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
        // calcolo lambda dell'algoritmo di gauss per le righe i e j
        float lambda = res[j][i] / res[i][i];

        for (int k = i; k < cols; k++) {
          res[j][k] -= lambda * res[i][k];
        }
      }
    }
    return res;
  }

  // funzione che calcola il rango della matrice
  int getRank() {
    // Iniziamo usando l'eliminazione di gauss
    Matrix<T> res = gauss();
    std::cout << "matrice res\n";

    // contiamo le righe non nulle
    int rank = 0;
    for (int i = 0; i < rows; i++) {
      bool null_row = false;
      for (int j = 0; j < cols; j++) {
        if (res[i][j] != 0) {
          null_row = true;
          break;
        }
      }
      if (null_row) {
        rank++;
      }
    }
    return rank;
  }

  T getDeterminant() {
    // Check if the matrix is square
    if (rows != cols) {
      // throw std::invalid_argument("matrice non quadrata");
      std::cout << "matrice non quadrata";
    }

    // uso gauss
    Matrix res = gauss();

    // il determinante è uguale al prodotto dei pivot
    T determinant = 1;
    for (int i = 0; i < rows; i++) {
      determinant *= res[i][i];
    }

    return determinant;
  }

  // calcola la trasposta della matrice (utile per calcolare l'inversa)
  Matrix transpose() {
    Matrix res(cols, rows);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        res[j][i] = matrix[i][j];
      }
    }
    return res;
  }

  // moltiplicazione riga per colonna
  Matrix mult(Matrix<T> B) {
    Matrix<T> res = Matrix(rows, B.getCols());

    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < B.getCols(); j++) {
        res[i][j] = 0;
        for (int k = 0; k < rows; k++) {
          res[i][j] += matrix[i][k] * B[k][j];
        }
      }
    }
    return res;
  }

  // override dell'operatore ostream
  friend std::ostream &operator<<(std::ostream &os, const Matrix m) {
    os << "{\n";
    for (size_t i = 0; i < m.getRows(); ++i) {
      os << " {";
      for (size_t j = 0; j < m.getCols(); ++j) {
        if (j > 0) {
          os << ", ";
        }
        os << m.matrix[i][j];
      }
      os << "}\n"; // non stampo nello standard output ma nello steam os
    }
    os << "}"; // non stampo nello standard output ma nello steam os

    return os; // restituisce riferimento puntatore
  }

  // ritorna la sottomatrice senza la riga r e la colonna c
  Matrix<T> submatrix(int r, int c) {
    // Sottomatrice
    Matrix<T> sub(rows - 1, cols - 1);

    // row ci indica a che riga siamo
    int row = 0;
    for (int i = 0; i < rows; i++) {
      // escludiamo la riga i
      if (i == r)
        continue;

      // col ci indica a che colonna siamo
      int col = 0;
      for (int j = 0; j < cols; j++) {
        // escludiamo la colonna c
        if (j == c)
          continue;

        // se non ci troviamo nella riga r o nella colonna c aggiungiamo
        // l'elemento
        sub[row][col] = matrix[i][j];
        col++;
      }
      row++;
    }

    return sub;
  }
  // override della moltiplicazione scalare
  Matrix<T> operator*(T b) {
    Matrix res(rows, cols);
    res.matrix = matrix;
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        res[i][j] *= b;
      }
    }
    return res;
  }

  // funzione che calcola l'inversa della matrice
  Matrix getInverse() {
    // se rows != cols non è invertibile
    if (rows != cols || rows != getRank()) {
      // errore
      std::cout << "errore, matrice non quadrata";
    }
    T det = getDeterminant();
    Matrix res(rows, cols);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < rows; j++) {
        Matrix sub = submatrix(i, j);
        res[i][j] = std::pow(-1, i + j) * sub.getDeterminant();
      }
    }
    // posso farlo perchè rows == cols
    res = res.transpose();
    T invdet = 1 / det;
    res = res * invdet;
    std::cout << "adesso printo res \n";
    std::cout << det;
    return res;
  }
};

int main() {
  srand(time(0)); // Seed per il random

  Matrix<int> m(4, 4);
  m[0][0] = -2;
  m[0][1] = -1;
  m[0][2] = 0;
  m[0][3] = 1;
  m[1][0] = -1;
  m[1][1] = 0;
  m[1][2] = 1;
  m[1][3] = 2;
  m[2][0] = 0;
  m[2][1] = 1;
  m[2][2] = 2;
  m[2][3] = 3;
  m[3][0] = 1;
  m[3][1] = 2;
  m[3][2] = 3;
  m[3][3] = 4;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      // se voglio inizializzare una matrice random
      // m[i][j] = rand() %10 +1;

      // inizializzazione di prova
      // m[i][j] = i + 1 + j;
    }
  }

  std::cout << "matrice m: \n" << m << std::endl;

  //std::cout << "matrice inversa: \n" << m.getInverse();

  std::cout << "matrice gauss: \n" << m.gauss();
  //std::cout << "determinante: " << m.getDeterminant() << "\n";
  //std::cout << "rango: " << m.getRank() << "\n";
  // std::cout << m.mult(m.getInverse());
  //std::cout << "matrice gauss: \n" << m.gauss();
  
}
