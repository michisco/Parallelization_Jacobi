#ifndef FFLOWJACOBI_H
#define FFLOWJACOBI_H

#include<stdlib.h>
#include<iostream>
#include<vector>
#include <ff/parallel_for.hpp>

#include "utimer.h"
#include "utilities.h"

using namespace std;

/**
 * @brief Base function that perform a parallel version of Jacobi algorithm with FastFlow.
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
 * @param time variable to store fastflow time
 * @return solution of Jacobi algorithm (last computation)
 */
vector<float> fflowJacobi(int maxIter, int matrixSize, int n_threads, vector<vector<float>> A, vector<float> b, long *time);

vector<float> fflowJacobi(int maxIter, int matrixSize, int n_threads, vector<vector<float>> A, vector<float> b, long *time){

    vector<float> old_value(matrixSize, 0);        //previous value of the computation
    vector<float> new_value(matrixSize, 0);        //new value of the computation

    //using a ParalleFor object with n_threads workers applying nonblocking policy
    ff::ParallelFor parallelCycle(n_threads, true);

    {
        utimer ff("Elapsed parallel_for time = ", time);

        for(int k=0; k<maxIter; k++){
            //apply parallelFor object
            parallelCycle.parallel_for(0, matrixSize, 1, 0, [&](const long i){
                float sum=0;
                for(int j=0; j<i; j++)
                    sum+=A[i][j]*old_value[j];
                for(int j=i+1; j<b.size(); j++)
                    sum+=A[i][j]*old_value[j];

                new_value[i]=(b[i]-sum)/A[i][i];
            }, n_threads);

            //check stopping criterion to stop the algorithm and save the number of iterations done
            if(checkStoppingCriteria(old_value, new_value))
                break;
            else
                old_value = new_value;
        }
    }

    return new_value;
}
#endif // FFLOWJACOBI_H
