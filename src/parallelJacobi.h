#ifndef PARALLELJACOBI_H
#define PARALLELJACOBI_H
#include<stdlib.h>
#include<iostream>
#include<vector>
#include <chrono>
#include <mutex>
#include <functional>
#include <atomic>
#include <thread>
#include <sched.h>

#include <barrier>

#include "utimer.h"
#include "utilities.h"

using namespace std;

/**
 * @brief Base function that perform a parallel version of Jacobi algorithm with barriers.
 *
 *        Get a square matrix and vector and compute the Jacobi method for determining
 *        the solution of a strictly diagonally dominant system of linear equation.
 *
 *        During the execution, calculate and store the time to perform the algorithm.
 *
 * @param maxIter maximum number of iterations
 * @param matrixSize dimension of matrix (nxn)
 * @param n_threads number of threads
 * @param A matrix
 * @param b vector
 * @param time variable to store parallel time
 * @return solution of Jacobi algorithm (last computation)
 */
vector<float> parallelJacobi(int maxIter, int matrixSize, int n_threads, vector<vector<float>> A, vector<float> b, long *time);

/**
 * @brief Base function that perform a parallel version of Jacobi algorithm with barriers using pinned threads.
 *
 *        Get a square matrix and vector and compute the Jacobi method for determining
 *        the solution of a strictly diagonally dominant system of linear equation.
 *
 *        During the execution, calculate and store the time to perform the algorithm.
 *
 * @param maxIter maximum number of iterations
 * @param matrixSize dimension of matrix (nxn)
 * @param n_threads number of threads
 * @param A matrix
 * @param b vector
 * @param time variable to store parallel time
 * @return solution of Jacobi algorithm (last computation)
 */
vector<float> parallelJacobiPinned(int maxIter, int matrixSize, int n_threads, vector<vector<float>> A, vector<float> b, long *time);

/**
 * @brief Base function that computes overhead of a parallel version of Jacobi algorithm with barrier.
 *
 * @param maxIter maximum number of iterations
 * @param matrixSize dimension of matrix (nxn)
 * @param n_threads number of threads
 * @return overhead time
 */
long computingOverhead(int maxIter, int matrixSize, int n_threads);

vector<float> parallelJacobi(int maxIter, int matrixSize, int n_threads, vector<vector<float>> A, vector<float> b, long *time){

    vector<float> old_value(matrixSize, 0);     //previous value of the computation
    vector<float> new_value(matrixSize, 0);     //new value of the computation

    int NrIter=maxIter;
    float sum_norm = 0;
    vector<float> norm(n_threads, 0);	//vector of the partial norms

    //barrier callback
    barrier barObj(n_threads, [&](){
    	//sum the partial norms
        for(int i=0; i < n_threads; i++){
            sum_norm += norm[i];
            norm[i]=0;
        }
        sum_norm=sum_norm/((float)(matrixSize));
        if(checkStoppingCriteria(sum_norm))
            NrIter=0;
        else{
            NrIter--;
            old_value=new_value;
            sum_norm=0;
        }
        return;
    });

    //thread lambda function
    auto sum=[&](int chunk_lower_bound, int chunk_upper_bound, int thread_i)	
    {
        while(NrIter>0){
            for (int i=chunk_lower_bound; i<chunk_upper_bound; i++){
                float sum=0;
                for(int j=0; j<i; j++)
                    sum += A[i][j]*old_value[j];
                for(int j=i+1; j<matrixSize; j++)
                    sum += A[i][j]*old_value[j];

                new_value[i]=(b[i]-sum)/A[i][i];

                //compute partial norm
                norm[thread_i]+=abs(old_value[i]-new_value[i]);
            }
            //call barrier
            barObj.arrive_and_wait();
        }
    };

    thread *t[n_threads];   //initialize n threads
    int n_chunk;    //chunks size
    n_chunk =  matrixSize % n_threads == 0  ? (matrixSize / n_threads) : (matrixSize / n_threads) + 1;

    int chunk_upper_bound = n_chunk;
    int chunk_lower_bound = 0;

    //main loop
    {
        utimer threadtime("Elapsed thread time = ", time);

        //assign sum lambda function with different chunk to each thread
        for(int thread_i=0; thread_i<n_threads; thread_i++){

            t[thread_i] = new thread(sum, chunk_lower_bound, min(chunk_upper_bound, matrixSize), thread_i);
            chunk_lower_bound = chunk_upper_bound;
            chunk_upper_bound = chunk_upper_bound+n_chunk;
        }

        for(int i=0; i<n_threads; i++)
            t[i]->join();
    }

    return new_value;
}

vector<float> parallelJacobiPinned(int maxIter, int matrixSize, int n_threads, vector<vector<float>> A, vector<float> b, long *time){

    vector<float> old_value(matrixSize, 0);     //previous value of the computation
    vector<float> new_value(matrixSize, 0);     //new value of the computation

    int NrIter=maxIter;
    float sum_norm = 0;
    vector<float> norm(n_threads, 0);	//vector of the partial norms

    //barrier lambda function
    barrier barObj(n_threads, [&](){
        for(int i=0; i < n_threads; i++){
            sum_norm += norm[i];
            norm[i]=0;
        }
        sum_norm=sum_norm/((float)(matrixSize));
        if(checkStoppingCriteria(sum_norm))
            NrIter=0;
        else{
            NrIter--;
            old_value=new_value;
            sum_norm=0;
        }
        return;
    });

    thread *t[n_threads];   //initialize n threads

    //thread lambda function
    auto sum=[&](int chunk_lower_bound, int chunk_upper_bound, int thread_i)
    {
    	//prepare data structure to store a set of cpus
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(thread_i, &cpuset);
        //sets the CPU affinity mask of a thread created 
        int rc = pthread_setaffinity_np(t[thread_i]->native_handle(), sizeof(cpu_set_t), &cpuset);

	//if returns a nonzero number then it is an error
        if (rc != 0)
            std::cerr << "Error calling pthread_setaffinity_np: " << rc << "\n";

        while(NrIter>0){
            for (int i=chunk_lower_bound; i<chunk_upper_bound; i++){
                float sum=0;
                for(int j=0; j<i; j++)
                    sum += A[i][j]*old_value[j];
                for(int j=i+1; j<matrixSize; j++)
                    sum += A[i][j]*old_value[j];

                new_value[i]=(b[i]-sum)/A[i][i];

                //compute partial norm
                norm[thread_i]+=abs(old_value[i]-new_value[i]);
            }
            //call barrier
            barObj.arrive_and_wait();
        }
    };

    int n_chunk;    //chunks size
    n_chunk =  matrixSize % n_threads == 0  ? (matrixSize / n_threads) : (matrixSize / n_threads) + 1;

    int chunk_upper_bound = n_chunk;
    int chunk_lower_bound = 0;

    //main loop
    {
        utimer pinnedthreadtime("Elapsed thread time = ", time);

        //assign sum lambda function with different chunk to each thread
        for(int thread_i=0; thread_i<n_threads; thread_i++){

            t[thread_i] = new thread(sum, chunk_lower_bound, min(chunk_upper_bound, matrixSize), thread_i);
            chunk_lower_bound = chunk_upper_bound;
            chunk_upper_bound = chunk_upper_bound+n_chunk;
        }

        for(int i=0; i<n_threads; i++)
            t[i]->join();
    }

    return new_value;
}

long computingOverhead(int maxIter, int matrixSize, int n_threads){
    long overhead_time;
    auto blank=[&](){return;};

    thread *t[n_threads];   //initialize n threads
    int n_chunk;    //chunks size
    n_chunk =  matrixSize % n_threads == 0  ? (matrixSize / n_threads) : (matrixSize / n_threads) + 1;

    int chunk_upper_bound = n_chunk;
    int chunk_lower_bound = 0;

    //main loop
    {
        utimer overheadtime("Elapsed overhead time = ", &overhead_time);

        //assign sum lambda function with different chunk to each thread
        for(int thread_i=0; thread_i<n_threads; thread_i++){

            t[thread_i] = new thread(blank);
            chunk_lower_bound = chunk_upper_bound;
            chunk_upper_bound = chunk_upper_bound+n_chunk;
        }

        for(int i=0; i<n_threads; i++)
            t[i]->join();
    }

    return overhead_time;
}

#endif // PARALLELJACOBI_H
