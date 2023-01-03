#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <ostream>
#include <vector>

//funzione per calcolare il massimo comune divisore con l'algoritmo di euclide
//un'alternativa è usare std::GCD(a,b) ma è supportata solo in C++17
int mcd(int a, int b) {
  if (b == 0)
    return a;
  return mcd(b, a % b);
}

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
    if (r == 0 || c == 0) {
      throw std::invalid_argument("matrice vuota");
    }
    // inizializzo la matrice (vuota) composta da r righe e c colonne. In
    // particolare ogni riga sarà un array di elementi
    matrix = std::vector<std::vector<T>>(r, std::vector<T>(c));
  }
  // implementiamo un costruttore "copia" che consente di creare un nuovo oggetto 
  // copiando gli attributi di quello passato come parametro
  Matrix(const Matrix &mat) {

    rows = mat.getRows();
    cols = mat.getCols();
    if (rows == 0 || cols == 0) {
      throw std::invalid_argument("matrice vuota");
    }
    matrix = mat.matrix;
  }

  // funzione che restituisce il numero di righe (la variabile è privata quindi
  // non è visualizzabile altrimenti)
  const int getRows() const { return rows; }

  // funzione che restituisce il numero di colonne (la variabile è privata
  // quindi non è visualizzabile altrimenti)
  const int getCols() const { return cols; }


  // funzione che applica alla matrice l'eliminazione di gauss

  const Matrix<T> gauss() {
    
    // usiamo una copia della matrice in modo tale da non modificare l'originale
    Matrix<T> res(*this);

    for (int i = 0; i < rows - 1; i++) {
      // flag serve a verificare la validità del pivot della riga i
      //se è true è tutto a posto
      bool flag = true;

      //se il pivot è nullo devo controllare che tutti gli elementi della colonna i sotto di lui siano nulli
      if (res[i][i] == 0) {
        flag = false;

        for (int j = i; j < rows; j++) {
          
          //se trovo una riga con un pivot valido, posso scambiare le righe i e j
          //a questo punto posso proseguire normalmente
          if (res[j][i] != 0) {
            flag = true;
            std::swap(res[j], res[i]);

            //siccome ora il pivot non è nullo posso uscire dal ciclo
            break;
          }
        }
      }
      // se res[i][i] ==0 e non ho trovato una riga da scambiare significa che
      // dalla riga i in giù, tutta la colonna i è vuota, quindi passo alla riga i+1
      if (!flag)
        continue;

      for (int j = i + 1; j < rows; j++) {
        // NOTA IMPORTANTE: per le matrici di interi calcoliamo comunque gauss,
        // ma non può essere utilizzato per il calcolo del determinante

        // nel caso delle matrici intere moltiplico la riga j per il pivot (per questo però non si potrà calcolare il determinante)
        // in modo tale da poter sottrarre alla riga j un multiplo della riga i garantendo di usare solo numeri interi
        if (std::is_integral<T>::value) {
          
          //k è l'indice della colonna
          for (int k = 0; k < cols; k++) {
            // questo algoritmo però spesso va in overflow anche senza usare numeri enormi, questo è inevitabile in alcuni casi
            // siccome c++ non gestisce gli overflow nelle operazioni tra interi, dobbiamo farlo manualmente
            T max =std::numeric_limits<T>::max()/res[i][i];
            if (max <0)
              max = -max;
            
            if (res[j][k]> max){
              //se vado in overflow alzo un'eccezione
              throw std::overflow_error("Errore: impossible utilizzare gauss con questa matrice, valori troppo grandi. Prova a convertire in un tipo più grande (ad esempio longint) oppure in float");
            }

            res[j][k] *= res[i][i];
          }
        }

        // calcolo lambda dell'algoritmo di gauss per le righe i e j
        float lambda = res[j][i] / res[i][i];

        //procediamo con l'eliminazione
        for (int k = i; k < cols; k++) {
          res[j][k] -= lambda * res[i][k];
        }

        //per le matrici di interi, semplifichiamo la riga j se possibile
        if (std::is_integral<T>::value) {

          //trovo il massimo comune divisore
          T max_div = res[j][0];
          for (int k = 1; k < cols; k++) {
            max_div = mcd(max_div, res[j][k]);
          }

          //nel caso di righe nulle max_div è 0, in quel caso lo poniamo =1
          if (max_div ==0)  
            max_div =1;
          //semplifico la riga
          for (int k = 1; k < cols; k++)
            res[j][k] /= max_div;

        }
      }
    }
    //ritorno la matrice res
    return res;
  }

  // funzione che calcola il rango della matrice
  const int getRank() {
    // Iniziamo usando l'eliminazione di gauss
    Matrix<T> res = Matrix(gauss());

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

  const T getDeterminant() {
    T det = 0;

    if (rows != cols) {
      throw std::invalid_argument("La matrice non è quadrata");
    } else {
      // caso base: la matrice è 1x1
      if (rows == 1) {
        return matrix[0][0];
      }
      // caso base: la matrice è 2x2
      if (rows == 2) {
        return matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];
      }
      // caso ricorsivo: calcola il determinante utilizzando il metodo di
      // Laplace
      else {
        for (int i = 0; i < cols; i++) {
          // calcolo il sottodeterminante utilizzando il metodo di Laplace sulla
          // sottomatrice eliminando la riga 0 e la colonna i
          Matrix<T> submatrix(rows - 1, cols - 1);
          for (int j = 1; j < rows; j++) {
            int submatrix_col = 0;
            for (int k = 0; k < cols; k++) {
              if (k == i)
                continue;
              submatrix[j - 1][submatrix_col] = matrix[j][k];
              submatrix_col++;
            }
          }
          T subdet = submatrix.getDeterminant();

          // aggiungo o sottraggo il sottodeterminante al risultato finale in
          // base al segno dell'elemento corrente (i, 0)
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
  const Matrix<T> mult(const Matrix<T> &B) {
    Matrix<T> res = Matrix(rows, B.cols);
    if (B.cols == 0 || B.rows == 0) {
      throw std::invalid_argument("Matrice B nulla");
    }
    if (cols != B.rows) {
      throw std::invalid_argument(
          "Errore: il numero di colonne della prima matrice deve essere uguale "
          "al numero di righe della seconda matrice");
    } else {
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
  const std::vector<T> mult_vect(const std::vector<T> &b) {
    std::vector<T> res(rows);

    if (cols != b.size()) {
      std::invalid_argument("Errore: il numero di colonne della matrice deve "
                            "essere uguale al numero di elementi del vettore");
    } else {
      for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
          res[i] += matrix[i][j] * b[j];
        }
      }
    }
    return res;
  }

  // Funzione per aggiungere un vettore colonna alla matrice
  const void AddColumn(std::vector<T> &column) {
    // Verifica che il vettore abbia lo stesso numero di elementi delle righe
    // della matrice
    if (column.size() != rows) {
      std::invalid_argument(
          "Errore: il vettore colonna ha un numero di elementi diverso dal "
          "numero di righe della matrice.");
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
    return column;
  }

  // risoluzione sistema lineare Ax=b    (dim x = num col A)

  std::vector<T> system(std::vector<T> &b) {
    if (std::is_integral<T>::value) {
      throw std::invalid_argument(
          "impossible risolvere il sistema lineare di interi");
    }
    std::vector<T> x(cols);

    // Verifica che il numero di elementi di b sia uguale al numero di righe
    // della matrice
    if (b.size() != rows) {
      throw std::invalid_argument(
          "Errore: il vettore b deve avere lo stesso numero di elementi del "
          "numero di righe della matrice");
    }

    else {
      // copio la matrice di partenza in orig
      Matrix<T> orig(rows, cols);
      orig = Matrix(*this);
      orig.AddColumn(b);
      Matrix<T> res = orig.gauss();
      Matrix<T> S(rows, cols); // copio res in S
      S = Matrix(res);

      std::vector<T> c = res.RemoveLastColumn(); // termini noti dopo gauss

      if (res.getDeterminant() == 0) {
        throw std::invalid_argument("le soluzioni non sono indipendenti");
      }

      else {
        // risolvo il sistema all'indietro
        for (int i = rows - 1; i >= 0; --i) {
          T sum = 0;
          for (int j = i + 1; j < cols; j++) {
            sum += S[i][j] * x[j];
          }
          x[i] = (c[i] - sum) / S[i][i];
        }
      }
    }
    return x;
  }

  Matrix<T> inverse() {
    if (std::is_integral<T>::value) {
      throw std::invalid_argument(
          "impossible calcolare l'inversa di una matrice di interi");
    }
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
      throw std::invalid_argument("matrice non invertibile");
    }

    else {

      Matrix<T> G = B.gauss();

      for (size_t i = 0; i < rows; ++i) {
        // Scegli il pivot

        for (size_t j = 0; j < i; ++j) {
          T factor = G[j][i] / G[i][i];
          for (size_t k = 0; k < 2 * cols; ++k) {
            G[j][k] -= G[i][k] * factor;
          }
        }
      }

      // Moltiplica ogni riga per il valore inverso del pivot
      for (size_t i = 0; i < rows; ++i) {
        T factor = 1.0 / G[i][i];
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

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &b) {
  os << "[";
  for (size_t i = 0; i < b.size(); ++i) {
    if (i > 0) {
      os << ",";
    }
    os << b[i];
  }
  os << "]";
  return os;
}

int main() {
  srand(time(0)); // Seed per il random
  int n = 15;
  Matrix<int> m(n, n);
  Matrix<int> m2(n, n);
  std::vector<int> v(n);
  for (int i = 0; i < n; i++) {
    v[i] = rand() % n + 1;
    for (int j = 0; j < n; j++) {
      // se voglio inizializzare una matrice random
      m[i][j] = rand() % 100 + 1;
      m2[i][j] = rand() % 100 + 1;

      // inizializzazione di prova
      // m[i][j] = i + 1 + j;
    }
  }
  m[0][0] = 0;
  std::cout << "matrice m: \n" << m << std::endl;
  std::cout << "vettore v: " << v << std::endl;
  try {
    std::cout << "matrice gauss: \n" << m.gauss();
    std::cout << "matrice inversa: \n" << m.inverse();
    // std::cout << "determinante: " << m.getDeterminant() << "\n";
    std::cout << "rango: " << m.getRank() << "\n";
    // std::cout << m.mult(m);
    std::cout << "matrice gauss: \n" << m.gauss();

    std::cout << "matrice m: \n" << m << std::endl;

    std::cout << m.system(v);

  } catch (std::exception &e) {
    std::cout << e.what();
  }
}
