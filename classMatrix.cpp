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

  /*Matrix(Matrix &mat) {
    rows = mat.getRows();
    cols = mat.getCols();
    matrix = mat.matrix;
  }*/

  // funzione che restituisce il numero di righe (la variabile è privata quindi
  // non è visualizzabile altrimenti)
  int getRows() const { return rows; }

  // funzione che restituisce il numero di colonne (la variabile è privata
  // quindi non è visualizzabile altrimenti)
  int getCols() const { return cols; }

  // funzione che restituisce la matrice dopo l'eliminazione di gauss
  // nota: la funzione restituisce una matrice float anche per le matrici di
  // interi
  Matrix<T> intGauss() {
    Matrix<T> res(rows, cols);
    res = *this;

    for (int i = 0; i < rows - 1; i++) {

      bool flag = true;
      std::cout << "res ii " << res[i][i] << " " << i << std::endl << res;
      if (res[i][i] == 0) {
        // flag per vedere se c'è una riga dove res[l][i] non è nullo
        flag = false;
        for (int l = i; l < rows; l++) {
          if (res[l][i] != 0) {
            flag = true;

            // scambio la riga l con la riga i
            std::cout << res;
            std::swap(res[l], res[i]);
            std::cout << res;
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
        T lambda = res[j][i] / res[i][i];

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
      throw std::invalid_argument("matrice non quadrata");
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

  // moltiplicazione riga - colonna matrice matrice
  Matrix<T> mult(const Matrix<T> &B) {
    Matrix<T> res = Matrix(rows, B.cols);

    if (cols != B.rows) {
      std::cout << "Errore: il numero di colonne della prima matrice deve "
                   "essere uguale al numero di righe della seconda matrice"
                << std::endl;
    } else {
      for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < B.cols; ++j) {
          for (size_t k = 0; k < rows; ++k) {
            // std::cout<< matrix[i][k] * B[k][j]<< "\n";
            res[i][j] += matrix[i][k] * B[k][j];
          }
        }
      }
      std::cout << "La matrice prodotto è: " << std::endl << res << std::endl;
    }
    return res;
  }

  // moltiplicazione matrice - vettore colonna
  std::vector<T> mult_vect(const std::vector<T> &b) {
    std::vector<T> res(rows);

    if (cols != b.size()) {
      std::cout << "Errore: il numero di colonne della matrice deve essere "
                   "uguale al numero di elementi del vettore"
                << std::endl;
    } else {
      for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
          res[i] += matrix[i][j] * b[j];
        }
      }
      std::cout << "Il vettore prodotto è: " << std::endl << res << std::endl;
    }
    return res;
  }

  // Funzione per aggiungere un vettore colonna alla matrice
  void AddColumn(std::vector<T> &column) {
    // Verifica che il vettore abbia lo stesso numero di elementi delle righe
    // della matrice
    if (column.size() != rows) {
      std::cout << "Errore: il vettore colonna ha un numero di elementi "
                   "diverso dal numero di righe della matrice."
                << std::endl;
    } else {
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
      matrix[i].pop_back(); // metodo pop_back per rimuovere l'ultimo elemento
                            // di ogni riga della matrice
    }
    cols = cols - 1;
    // std::cout << column << std::endl;
    // std::cout << matrix;
    return column;
  }

  // risoluzione sistema lineare Ax=b    (dim x = num col A)

  // funzione per l'eliminazione di Gauss
  std::vector<T> system(std::vector<T> &b) {
    std::vector<T> x(cols);

    // Verifica che il numero di elementi di b sia uguale al numero di righe
    // della matrice
    if (b.size() != rows) {
      // throw std::invalid_argument("Errore: il vettore b deve avere lo stesso
      // numero di elementi del numero di righe della matrice");
      std::cout << "Errore: il vettore b deve avere lo stesso numero di "
                   "elementi del numero di righe della matrice"
                << std::endl;
    }

    else {
      // copio la matrice di partenza in orig
      Matrix<T> orig(rows, cols);
      for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
          orig[i][j] = matrix[i][j];
        }
      }
      orig.AddColumn(b);
      Matrix<T> res = orig.gauss();
      Matrix<T> S(rows, cols); // copio res in S
      for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
          S[i][j] = res[i][j];
        }
      }

      std::vector<T> c = res.RemoveLastColumn(); // termini noti dopo gauss

      if (res.getDeterminant() == 0) {
        std::cout << "Stop: le soluzioni non sono indipendenti" << std::endl;
      }

      else {
        std::cout << "Avanti: " << std::endl;
        // risolvo il sistema all'indietro
        for (int i = rows - 1; i >= 0; --i) {
          T sum = 0;
          for (int j = i + 1; j < cols; j++) {
            sum += S[i][j] * x[j];
          }
          x[i] = (c[i] - sum) / S[i][i];
        }

        // stampo il vettore delle soluzioni
        // std::cout << "Il vettore delle soluzioni è: " << x << std::endl;
        for (int i = 0; i < cols; ++i) {
          std::cout << "Il vettore delle soluzioni è: " << x[i] << std::endl;
        }
      }
    }
    return x;
  }

  Matrix<T> inverse() {
    Matrix<T> inv(rows, cols);

    // matrirce identità
    Matrix<T> id(rows, cols);
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        if (i == j) {
          id[i][j] = 1;
        } else {
          id[i][j] = 0;
        }
      }
    }
    // id.print();  std::cout << std::endl;

    Matrix<T> B(rows, 2 * cols);
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        B[i][j] = matrix[i][j];
      }
      for (int j = cols; j < 2 * cols; j++) {
        B[i][j] = id[i][j - cols];
      }
    }

    if (rows != cols) { // | A.getDeterminant() == 0) {
      // throw std::invalid_argument("matrice non invertibile");
      std::cout << "Errore: la matrice non è invertibile" << std::endl;
    }

    else {
      std::cout<<B;
      std::cout << std::endl;
      Matrix<T> G = B.gauss();
      std::cout << std::endl;

      for (size_t i = 0; i < rows; ++i) {
        // Scegli il pivot

        for (size_t j = 0; j < i; ++j) {
          double factor = G[j][i] / G[i][i];
          for (size_t k = 0; k < 2 * cols; ++k) {
            G[j][k] -= G[i][k] * factor;
          }
        }
      }

      // Moltiplica ogni riga per il valore inverso del pivot
      for (size_t i = 0; i < rows; ++i) {
        double factor = 1.0 / G[i][i];
        for (size_t j = 0; j < 2 * rows; j++) {
          G[i][j] *= factor;
        }
      }

      // Estrai la matrice inversa dalla destra della matrice estesa
      for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
          inv[i][j] = G[i][j + cols];
        }
      }
      std::cout << inv;
    }

    return inv;
  }

  // override dell'operatore [][] che permette di scrivere Matrix[i][j]
  inline std::vector<T> &operator[](size_t i) { return matrix[i]; }

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



int main() {
  srand(time(0)); // Seed per il random

  Matrix<float> m(4, 4);
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
  std::cout << "matrice inversa: \n" << m.inverse();
  std::cout << "matrice gauss: \n" << m.gauss();

  std::cout << "determinante: " << m.getDeterminant() << "\n";
  std::cout << "rango: " << m.getRank() << "\n";
  std::cout << m.mult(m.inverse());
  std::cout << "matrice gauss: \n" << m.gauss();
}
