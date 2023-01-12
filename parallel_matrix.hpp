#ifndef __parallel_matrix_hpp__  
#define __parallel_matrix_hpp__ 

#include <iostream>
#include <ostream>
#include <exception>

#include "opencl_processor.hpp"

#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120

#include <CL/opencl.hpp> //per binding C++ di OpenCL

/**
 * @brief Una classe per la rappresentazione delle matrici tramite programmazione parallela
*/
class ParallelMatrix : public OpenCLProcessor {
protected:
    std::string kernel_name;   //!< nome della funzione del kernel

    /**
     * @brief Accoda il kernel 
     * 
     * @param queue è una coda di comandi
     * @param kernel è il kernel da accodare
     * @param work_item_size è il numero dei work-item da attivare
     */
    virtual void enqueue_kernel(cl::CommandQueue& queue, cl::Kernel& kernel, 
                                const cl::NDRange& work_item_size) {
        // sottometto il kernel
        queue.enqueueNDRangeKernel(kernel, cl::NullRange, work_item_size);
    }

public:
    /**
     * @brief Costruttore
     */
    ParallelMatrix(): OpenCLProcessor()
    {   
        // carico il kernel
        auto kernel_source = load_file("file.cl");
        // compilo il programma
	    build_program(kernel_source);
    }

    /**
     * @brief Nome della classe
     * 
     * @return il nome della classe
     */
    static std::string name()
    {
        return "ParallelMatrix";
    }

    /**
     * @brief Calcola il rango di una matrice
     * 
     * @param A la matrice di partenza
     * @param rows è il numero di righe di A
     * @param cols è il numero di colonne di A
     * @param rank è il rango di A
     */
    void getRank(float* A, const unsigned int rows, const unsigned int cols, int* rank)
    {
        // preparo i buffer
        const unsigned int datasize_A = rows*cols*sizeof(float);
        const unsigned int datasize_rank = sizeof(int);
        
        const auto input_flags = CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS | CL_MEM_COPY_HOST_PTR;
        const auto output_flags = CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR;

        cl::Buffer aBuf(context, input_flags, datasize_A);
        cl::Buffer rankBuf(context, output_flags, datasize_rank);

        // preparo il kernel e i suo parametri attuali
        cl::Kernel kernel(program, "getRank");
        kernel.setArg(0, aBuf);
        kernel.setArg(1, rows);
        kernel.setArg(2, cols);
        kernel.setArg(3, rankBuf);

        // creo una coda di comandi
        cl::CommandQueue queue(context, device);

        // sottometto il kernel
        enqueue_kernel(queue, kernel, cl::NDRange(rows,cols));

        queue.enqueueReadBuffer(rankBuf, CL_TRUE, 0, datasize_rank, rank);
    }

     
    /**
     * @brief Calcola il determinante di una matrice
     *
     * @param A la matrice di partenza
     * @param rows è il numero di righe di A 
     * @param cols è il numero di colonne di A 
     * @param determinant è il determinante di A 
     */
    void getDeterminant(float* A, const unsigned int rows, 
                     const unsigned int cols, float* determinant) 
    {
        if (rows != cols) {
            throw std::invalid_argument("La matrice non è quadrata!");
        }

        // preparo i buffer
        const unsigned int datasize_A = rows*cols*sizeof(float);
        const unsigned int datasize_determinant = sizeof(float);
        
        const auto input_flags = CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS | CL_MEM_COPY_HOST_PTR;
        const auto output_flags = CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR;

        cl::Buffer aBuf(context, input_flags, datasize_A);
        cl::Buffer detBuf(context, output_flags, datasize_determinant);

        // preparo il kernel e i suo parametri attuali
        cl::Kernel kernel(program, "getDeterminant");
        kernel.setArg(0, aBuf);
        kernel.setArg(1, rows);
        kernel.setArg(2, cols);
        kernel.setArg(3, detBuf);

        // creo una coda di comandi
        cl::CommandQueue queue(context, device);

        // sottometto il kernel
        enqueue_kernel(queue, kernel, cl::NDRange(rows,cols));

        queue.enqueueReadBuffer(detBuf, CL_TRUE, 0, datasize_determinant, determinant);
    }


    /**
     * @brief Calcola il prodotto tra due matrici
     * 
     * @param A è la prima matrice
     * @param B è la seconda matrice
     * @param C è la matrice risultato
     * @param A_rows è il numero di righe di A
     * @param A_cols è il numero di colonne di A e il numero di righe di B
     * @param B_cols è il numero di colonne di B
     */
    void mult_matrix(const float* A, const float* B,  float* C,
		             const unsigned int A_rows, const unsigned int A_cols, const unsigned int B_cols)
    {   
        // preparo i buffer
        const unsigned int datasize_A = A_rows*A_cols*sizeof(float);
        const unsigned int datasize_B = A_cols*B_cols*sizeof(float);
        const unsigned int datasize_C = A_rows*B_cols*sizeof(float);
        
        const auto input_flags = CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS | CL_MEM_COPY_HOST_PTR;
        const auto output_flags = CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR;

        cl::Buffer aBuf(context, input_flags, datasize_A, const_cast<float *>(A));
        cl::Buffer bBuf(context, input_flags, datasize_B, const_cast<float *>(B));
        cl::Buffer cBuf(context, output_flags, datasize_C);

        // preparo il kernel e i suo parametri attuali
        cl::Kernel kernel(program,"mult_matrix");
        kernel.setArg(0, aBuf);
        kernel.setArg(1, bBuf);
        kernel.setArg(2, cBuf);
        kernel.setArg(3, A_rows);
        kernel.setArg(4, A_cols);
        kernel.setArg(5, B_cols);

        // creo una coda di comandi
        cl::CommandQueue queue(context, device);

        // sottometto il kernel
        enqueue_kernel(queue, kernel, cl::NDRange(A_rows,B_cols));

        queue.enqueueReadBuffer(cBuf, CL_TRUE, 0, datasize_C, C);
    }


    /**
     * @brief Calcola il prodotto tra una matrice e un vettore
     * 
     * @param A è la matrice
     * @param b è il vettore
     * @param res è il vettore risultato
     * @param rows è il numero di righe di A
     * @param cols è il numero di colonne di A e la dimensione di b
     */
    void mult_vect(const float* A, const float* b,  float* res,
		 const unsigned int rows, const unsigned int cols)
    {   
        // preparo i buffer
        const unsigned int datasize_A = rows*cols*sizeof(float);
        const unsigned int datasize_b = cols*sizeof(float);
        const unsigned int datasize_res = rows*sizeof(float);
        
        const auto input_flags = CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS | CL_MEM_COPY_HOST_PTR;
        const auto output_flags = CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR;

        cl::Buffer aBuf(context, input_flags, datasize_A, const_cast<float *>(A));
        cl::Buffer bBuf(context, input_flags, datasize_b, const_cast<float *>(b));
        cl::Buffer resBuf(context, output_flags, datasize_res);

        // preparo il kernel e i suo parametri attuali
        cl::Kernel kernel(program, "mult_vect");
        kernel.setArg(0, aBuf);
        kernel.setArg(1, bBuf);
        kernel.setArg(2, resBuf);
        kernel.setArg(3, rows);
        kernel.setArg(4, cols);

        // creo una coda di comandi
        cl::CommandQueue queue(context, device);

        // sottometto il kernel
        enqueue_kernel(queue, kernel, cl::NDRange(rows,cols));     

        queue.enqueueReadBuffer(resBuf, CL_TRUE, 0, datasize_res, res);
    }


    /**
     * @brief Risolve il sistema lineare Ax=b con il metodo di Gauss-Jordan
     * 
     * @param A è la matrice
     * @param b è il vettore dei termini noti
     * @param x è il vettore risultato
     * @param rows è il numero di righe di A e la dimensione di b
     * @param cols è il numero di colonne di A
     */
    void system(float* A, float* b, float* x, unsigned int rows, unsigned int cols)
    {
        if (rows!=cols) {
            throw std::invalid_argument("Stop: le soluzioni non sono indipendenti");
        }

        // preparo i buffer
        const unsigned int datasize_A = rows*cols*sizeof(float);
        const unsigned int datasize_b = rows*sizeof(float);
        const unsigned int datasize_x = cols*sizeof(float);
        
        const auto input_flags = CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS | CL_MEM_COPY_HOST_PTR;
        const auto output_flags = CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR;

        cl::Buffer aBuf(context, input_flags, datasize_A);
        cl::Buffer bBuf(context, input_flags, datasize_b);
        cl::Buffer xBuf(context, output_flags, datasize_x);

        // preparo il kernel e i suo parametri attuali
        cl::Kernel kernel(program, "system");
        kernel.setArg(0, aBuf);
        kernel.setArg(1, bBuf);
        kernel.setArg(2, xBuf);
        kernel.setArg(3, rows);
        kernel.setArg(4, cols);

        // creo una coda di comandi
        cl::CommandQueue queue(context, device);

        // sottometto il kernel
        enqueue_kernel(queue, kernel, cl::NDRange(rows,cols));

        queue.enqueueReadBuffer(xBuf, CL_TRUE, 0, datasize_x, x);
    }


    /**
     * @brief Calcola l'inversa di una matrice 
     * 
     * @param A è la matrice
     * @param inv è la matrice inversa
     * @param rows è il numero di righe di A
     * @param cols è il numero di colonne di A
     */
    void  inverse(float* A, float* inv, unsigned int rows, unsigned int cols)
    {
        if (rows!=cols) {
            throw std::invalid_argument("Errore: la matrice non è invertibile");
        }

        // preparo i buffer
        const unsigned int datasize_A = rows*cols*sizeof(float);
        const unsigned int datasize_inv = rows*cols*sizeof(float);
        
        const auto input_flags = CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS | CL_MEM_COPY_HOST_PTR;
        const auto output_flags = CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY | CL_MEM_COPY_HOST_PTR;

        cl::Buffer aBuf(context, input_flags, datasize_A);
        cl::Buffer invBuf(context, output_flags, datasize_inv);

        // preparo il kernel e i suo parametri attuali
        cl::Kernel kernel(program, "inverse");
        kernel.setArg(0, aBuf);
        kernel.setArg(1, invBuf);
        kernel.setArg(2, rows);
        kernel.setArg(3, cols);

        // creo una coda di comandi
        cl::CommandQueue queue(context, device);

        // sottometto il kernel
        enqueue_kernel(queue, kernel, cl::NDRange(rows,cols));

        queue.enqueueReadBuffer(invBuf, CL_TRUE, 0, datasize_inv, inv);
    }
};


#endif  // __parallel_matrix_hpp__
