#ifndef SEQUENTIALJACOBI_H
#define SEQUENTIALJACOBI_H
#include<stdlib.h>
#include<iostream>
#include<vector>

#include "utimer.h"
#include "utilities.h"

using namespace std;

/**
 * @brief Base function that perform a sequential version of Jacobi algorithm.
 *
 *        Get a square matrix and vector and compute the Jacobi method for determining
 *        the solution of a strictly diagonally dominant system of linear equation.
 *
 *        During the execution, calculate and store the time to perform the algorithm
 *        and the number of iterations done.
 *
 * @param maxIter maximum number of iterations
 * @param matrixSize dimesion of matrix (nxn)
 * @param A matrix
 * @param b vector
 * @param time variable to store sequential time
 * @param nrIter number of iterations done
 * @return solution of Jacobi algorithm (last computation)
 */
vector<float> seqJacobi(int maxIter, int matrixSize, vector<vector<float>> A, vector<float> b, long *time, int *nrIter);

vector<float> seqJacobi(int maxIter, int matrixSize, vector<vector<float>> A, vector<float> b, long *time, int *nrIter){

    vector<float> old_value(matrixSize, 0);        //previous value of the computation
    vector<float> new_value(matrixSize, 0);        //new value of the computation

    float sum;
    utimer seq("Elapsed sequencial time = ", time);

    //iterative Jacobi algorithm
    for(int iter=0; iter<maxIter; iter++){
        for(int i = 0; i < matrixSize; i++){
            sum = 0;
            for(int j = 0; j < i; j++)
            	sum += A[i][j]*old_value[j];
            for(int j = i+1; j < matrixSize; j++)
            	sum += A[i][j]*old_value[j];

            new_value[i] = (b[i] - sum) / A[i][i];
        }

        //check stopping criterion to stop the algorithm and save the number of iterations done
        if(checkStoppingCriteria(old_value, new_value)){
            *nrIter = iter;
            break;
        }
        else
            old_value = new_value;
    }

    return new_value;
}
#endif // SEQUENTIALJACOBI_H
