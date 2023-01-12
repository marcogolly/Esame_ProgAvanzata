# Esame di Programmazione Avanzata e Parrallela
Cancian Piero, 
Gollinucci Marco, 
Majer William, 
Vattolo Carlotta

## Il Progetto
Il progetto riguarda l'implementazione di una classe [Matrix](/matrix.hpp) in linguaggio C++ per la rappresentazione delle matrici. All'interno della classe stessa vengono definiti i seguenti metodi:
* *getRank* per il calcolo del rango della matrice utilizzando il metodo di Gauss
* *getDeterminant* per il calcolo del determinante utilizzando il metodo di Laplace
* *mult_matrix* per calcolare la moltiplocazione righe-colonne di due matrici
* *mult_vect* per calcolare la moltiplicazione tra una matrice ed un vettore
* *system* per risolvere un sistema lineare del tipo Ax=b attraverso l'eliminazione di Gauss
* *inverse* per lacolcolare l'inversa di una matrice con il metodo di Gauss-Jordan

Oltre a questi metodi principali, vengono implementati altri metodi ausiliari quali:
* *getRows* per ottenere il numero di righe
* *getCols* per ottenere il numero di colonne
* *gauss* per calcolare la matrice ridotta con l'eliminazione di Gauss
* *RemoveLastColumn* per rimuovere l'ultima colonna
* *AddColumn* per aggiungere un vettore colonna
* *print* per la stampa della matrice

Esterni alla classe vengono poi implemetate delle funzioni e degli operatri utili al fine di snellire il codice.
Viene poi implementata la medesima classe ma rinominandola [ParallelMatrix](/parallel_matrix.hpp) utilizzando OpenCL per confrontare i tempi di esecuzione del codice.

## I file

elenco file con breve descrizione

requisiti, versione oepncl verisone c++

come si compila (comandi)

come si runna il main

piccol esposizione dei risultati sui garfici

c++ da 10x10 in poi molto lento con det, system, inverse
