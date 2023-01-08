//#ifndef __matrici_hpp__  
//#define __matrici_hpp__ 

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <ostream>
#include <vector>
#include "opencl_processor.hpp"

#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120

#include <CL/opencl.hpp> //per binding C++ di OpenCL


/*-------------------------------------- VERSIONE DI RISERVA -------------------------------------
//scrivo la classe per selezionare il dispositivo sul quale vogliamo
//mandare in esecuzione qaulcosa e compilare il programma
class OpenCLProcessor{
protected:    
    cl::Program program;
    cl::Device device;
    cl::Contex context; 
public:
    OpenCLProcessor(){
        //dovrà trovare tutte le piattaforme disponibili sulla macchina 
        //in cui sto cercando di istanziare un oggetto di questo tipo
            std::vector<cl::Platform> platforms;//vettore delle piattaforme 
            cl::Platform::get(&platforms);//per interrogare il sistema 
            //per sapere quali sono le piattaforme disponibili
            if(platforms.empty()){
                //se non ci sono piattaforme lancio un eccezione
                throw std::runtime_error("Nessuna piattaforma OpenCL disponibile!");
            }
            //se il vettore contiene almeno un elemento allora quel elemento 
            //è una piattaforma disponibile sulla mia macchina
            //mi chieod quali sono i dispositivi disponibili su quella piattaforma
            std::vector<cl::Device> devices;
            platforms[0].getDevices(CL_DEVICE_TYPE_ALL, &devices);
            if(devices.empty()){
                //se non ci sono device lancio un eccezione
                throw std::runtime_error("Nessun dispositivo OpenCL disponibile nella piattaforma!");
            }
            device = devices[0];
            context = cl::Context(device); 
    }   
    //costruisco un metodo che ci consente di caricare un programma all'interno del membro
    void build_program(const std::string kernel_source){
        //carico il sorgente
        cl::Program::Sources sources{kernel_source};
        //carico il programma
        program= cl::Program(context, sources);
        if (program.build() != CL_BUILD_SUCCES){
            std::ostringstream oss;
            oss<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
            throw std::runtime_error(oss.str());
        }
    }
};
//scrivo una funzione che servirà per caricare da un file il sorgente del kernell
std::string load_content(const std::string filename)
{
    //classe per leggere i file di output
    std::ifstream file(filename);
    std::string constent(std::istreambuf_iterator<char>(file), //iteratore inizio file
                        (std::istreambuf_iterator<char>()));
    return content;                    
} 
//scrivo la classe che caricherà il kernel e i suoi parametri
class ParallelMatrix : public OpenCLProcessor{
public:
    //costruttore classe
    ParallelMatrixMult() : OpenCLProcessor()
    {
        auto kernel_source = load_content("file.cl");
        build_program(kernel_source);
    }
        
}*/


/*-------------------------------- VERSIONE TIPO PROF --------------------------------*/

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
     * @brief Costruttori
     */
    
    ParallelMatrix(const std::string& kernel_name1)
        : OpenCLProcessor(), kernel_name("rank")
    {
        auto kernel_source = load_file("file.cl");
	    build_program(kernel_source);
    }

    ParallelMatrix(const std::string& kernel_name2)
        : OpenCLProcessor(), kernel_name("determinant")
    {
        auto kernel_source = load_file("file.cl");
	    build_program(kernel_source);
    }

    ParallelMatrix(const std::string& kernel_name3)
        : OpenCLProcessor(), kernel_name("mult_matrix")
    {
        auto kernel_source = load_file("file.cl");
	    build_program(kernel_source);
    }

    ParallelMatrix(const std::string& kernel_name4)
        : OpenCLProcessor(), kernel_name("mult_vect")
    {
        auto kernel_source = load_file("file.cl");
	    build_program(kernel_source);
    }

    ParallelMatrix(const std::string& kernel_name5)
        : OpenCLProcessor(), kernel_name("system")
    {
        auto kernel_source = load_file("file.cl");
	    build_program(kernel_source);
    }

    ParallelMatrix(const std::string& kernel_name6)
        : OpenCLProcessor(), kernel_name("inverse")
    {
        auto kernel_source = load_file("file.cl");
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

    void rank(const float* A, const unsigned int rows, const unsigned int cols, int* rank)
    {
        // preparo i buffer
        const unsigned int datasize_A = rows*cols*sizeof(float);
        const unsigned int datasize_rank = sizeof(int);
        
        const auto input_flags = CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS | CL_MEM_COPY_HOST_PTR;
        const auto output_flags = CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY;

        cl::Buffer aBuf(context, input_flags, datasize_A, const_cast<float *>(A));
        cl::Buffer rankBuf(context, input_flags, datasize_rank);

        // preparo il kernel e i suo parametri attuali
        cl::Kernel kernel(program, kernel_name.c_str());
        kernel.setArg(0, aBuf);
        kernel.setArg(1, rows);
        kernel.setArg(2, cols);
        kernel.setArg(3, rankBuf);

        // creo una coda di comandi
        cl::CommandQueue queue(context, device);

        // sottometto il kernel
        enqueue_kernel(queue, kernel, cl::NDRange(sizeof(int)));

        // A
        queue.enqueueReadBuffer(rankBuf, CL_TRUE, 0, datasize_rank, rank);
    }

    void determinant(const float* A, const unsigned int rows, 
                      const unsigned int cols, float* determinant) 
    {
        // preparo i buffer
        const unsigned int datasize_A = rows*cols*sizeof(float);
        const unsigned int datasize_determinant = sizeof(float);
        
        const auto input_flags = CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS | CL_MEM_COPY_HOST_PTR;
        const auto output_flags = CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY;

        cl::Buffer aBuf(context, input_flags, datasize_A, const_cast<float *>(A));
        cl::Buffer determinantBuf(context, input_flags, datasize_determinant);

        // preparo il kernel e i suo parametri attuali
        cl::Kernel kernel(program, kernel_name.c_str());
        kernel.setArg(0, aBuf);
        kernel.setArg(1, rows);
        kernel.setArg(2, cols);
        kernel.setArg(3, determinantBuf);

        // creo una coda di comandi
        cl::CommandQueue queue(context, device);

        // sottometto il kernel
        enqueue_kernel(queue, kernel, cl::NDRange(sizeof(int)));

        // A
        queue.enqueueReadBuffer(determinantBuf, CL_TRUE, 0, datasize_determinant, determinant);
    }

    void mult_matrix(const float* A, const float* B,  float* C,
		             const unsigned int A_rows, const unsigned int A_cols, const unsigned int B_cols)
    {
        // preparo i buffer
        const unsigned int datasize_A = A_rows*A_cols*sizeof(float);
        const unsigned int datasize_B = A_cols*B_cols*sizeof(float);
        const unsigned int datasize_C = A_rows*B_cols*sizeof(float);
        
        const auto input_flags = CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS | CL_MEM_COPY_HOST_PTR;
        const auto output_flags = CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY;

        cl::Buffer aBuf(context, input_flags, datasize_A, const_cast<float *>(A));
        cl::Buffer bBuf(context, input_flags, datasize_B, const_cast<float *>(B));
        cl::Buffer cBuf(context, input_flags, datasize_C);

        // preparo il kernel e i suo parametri attuali
        cl::Kernel kernel(program, kernel_name.c_str());
        kernel.setArg(0, aBuf);
        kernel.setArg(1, bBuf);
        kernel.setArg(2, cBuf);
        kernel.setArg(3, A_rows);
        kernel.setArg(4, A_cols);
        kernel.setArg(5, B_cols);

        // creo una coda di comandi
        cl::CommandQueue queue(context, device);

        // sottometto il kernel
        enqueue_kernel(queue, kernel, cl::NDRange(A_rows*B_cols));

        // A
        queue.enqueueReadBuffer(cBuf, CL_TRUE, 0, datasize_C, C);
    }

    /**
     * @param A è la prima matrice
     * @param B è il vettore
     * @param rows è il numero di righe di A
     * @param cols è il numero di colonne di A
     * @param M è la dimensione di B
     */
    void mult_vect(const float* A, const float* b,  float* res,
		 const unsigned int rows, const unsigned int cols)
    {
        // preparo i buffer
        const unsigned int datasize_A = rows*cols*sizeof(float);
        const unsigned int datasize_b = cols*sizeof(float);
        const unsigned int datasize_res = rows*sizeof(float);
        
        const auto input_flags = CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS | CL_MEM_COPY_HOST_PTR;
        const auto output_flags = CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY;

        cl::Buffer aBuf(context, input_flags, datasize_A, const_cast<float *>(A));
        cl::Buffer bBuf(context, input_flags, datasize_b, const_cast<float *>(b));
        cl::Buffer resBuf(context, input_flags, datasize_res);

        // preparo il kernel e i suo parametri attuali
        cl::Kernel kernel(program, kernel_name.c_str());
        kernel.setArg(0, aBuf);
        kernel.setArg(1, bBuf);
        kernel.setArg(2, resBuf);
        kernel.setArg(3, rows);
        kernel.setArg(4, cols);

        // creo una coda di comandi
        cl::CommandQueue queue(context, device);

        // sottometto il kernel
        enqueue_kernel(queue, kernel, cl::NDRange(rows));

        // A
        queue.enqueueReadBuffer(resBuf, CL_TRUE, 0, datasize_res, res);
    }

    void system(float* A, float* b, float* x, unsigned int rows, unsigned int cols)
    {
        // preparo i buffer
        const unsigned int datasize_A = rows*cols*sizeof(float);
        const unsigned int datasize_b = rows*sizeof(float);
        const unsigned int datasize_x = cols*sizeof(float);
        
        const auto input_flags = CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS | CL_MEM_COPY_HOST_PTR;
        const auto output_flags = CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY;

        cl::Buffer aBuf(context, input_flags, datasize_A, const_cast<float *>(A));
        cl::Buffer bBuf(context, input_flags, datasize_b, const_cast<float *>(b));
        cl::Buffer xBuf(context, input_flags, datasize_x);

        // preparo il kernel e i suo parametri attuali
        cl::Kernel kernel(program, kernel_name.c_str());
        kernel.setArg(0, aBuf);
        kernel.setArg(1, bBuf);
        kernel.setArg(2, xBuf);
        kernel.setArg(3, rows);
        kernel.setArg(4, cols);

        // creo una coda di comandi
        cl::CommandQueue queue(context, device);

        // sottometto il kernel
        enqueue_kernel(queue, kernel, cl::NDRange(cols));

        // A
        queue.enqueueReadBuffer(xBuf, CL_TRUE, 0, datasize_x, x);
    }


    void  inverse(float* A, float* inv, unsigned int N)
    {
        // preparo i buffer
        const unsigned int datasize_A = N*N*sizeof(float);
        const unsigned int datasize_inv = N*N*sizeof(float);
        
        const auto input_flags = CL_MEM_READ_ONLY | CL_MEM_HOST_NO_ACCESS | CL_MEM_COPY_HOST_PTR;
        const auto output_flags = CL_MEM_READ_WRITE | CL_MEM_HOST_READ_ONLY;

        cl::Buffer aBuf(context, input_flags, datasize_A, const_cast<float *>(A));
        cl::Buffer invBuf(context, input_flags, datasize_inv);

        // preparo il kernel e i suo parametri attuali
        cl::Kernel kernel(program, kernel_name.c_str());
        kernel.setArg(0, aBuf);
        kernel.setArg(1, invBuf);
        kernel.setArg(2, N);

        // creo una coda di comandi
        cl::CommandQueue queue(context, device);

        // sottometto il kernel
        enqueue_kernel(queue, kernel, cl::NDRange(N*N));

        // A
        queue.enqueueReadBuffer(invBuf, CL_TRUE, 0, datasize_inv, inv);
    }
};



