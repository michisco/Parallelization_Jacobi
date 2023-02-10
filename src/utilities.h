#ifndef UTILITIES_H
#define UTILITIES_H
#pragma once

#include <cmath>
#include <stdlib.h>
#include <vector>
#include <iostream>

using namespace std;

double EPSILON = 1e-11;
float MIN_VALUE = 0;
float MAX_VALUE = 3;

/**
 * @brief Check if the two vectors reached the stopping criterion computing the norm.
 *
 * @param a old value vector
 * @param b new value vector
 * @return a boolean value that shows if the stopping criterion occurs
 */
bool checkStoppingCriteria(vector<float> a, vector<float> b);

/**
 * @brief Compute the speedup
 * @param tseq sequential time
 * @param tpar parallel time
 * @return a float value that represents the speedup of algorithm
 */
float speedup(long tseq, long tpar);

/**
 * @brief Compute the efficiency
 * @param tseq sequential time
 * @param tpar parallel time
 * @param n_threads numbero of threads
 * @return a float value that represents the efficiency of algorithm
 */
float efficiency(long tseq, long tpar, int n_threads);

/**
 * @brief Compute the scalability
 * @param tpar1 parallel time with 1 thread
 * @param tparn parallel time with n thread
 * @return a float value that represents the scalability of algorithm
 */
float scalability(long tpar1, long tparn);

/**
* @brief Generate a random sizexsize matrix
* @param size dimension of matrix (nxn)
* @return a float random square matrix sizexsize.
*/
vector<vector<float>> matrixGenerator(int size);

/**
* @brief Generate a random right side vector
* @param size dimension of vector (n)
* @return a float right side vector with dimension 'size'
*/
vector<float> RHSVectorGenerator(int size);

/**
* @brief Print a matrix
* @param M float matrix
*/
void printMatrix(vector<vector<float>> M);

/**
* @brief Print a vector, the result of Jacobi algorithm
* @param b float vector
*/
void printResult(vector<float> b);

bool checkStoppingCriteria(vector<float> a, vector<float> b){
    float norm = 0;
    for(int i = 0; i < a.size(); i++)
        norm += abs(a[i] - b[i]);
    norm = norm / (float) a.size();

    return (norm <= EPSILON);
}

float speedup(long tseq, long tpar){
    return ((float)(tseq)/tpar);
}

float efficiency(long tseq, long tpar, int n_threads){
    return ((float)(tseq)/(n_threads * tpar));
}

float scalability(long tpar1, long tparn){
    return ((float)(tpar1)/tparn);
}

/**
 * @brief Check if the algorithm reached the stopping criterion using a pre-computed norm.
 *
 * @param norm matrix norm
 * @return a boolean value that shows if the stopping criterion occurs
 */
bool checkStoppingCriteria(float norm){
    return (norm <= EPSILON);
}

vector<vector<float>> matrixGenerator(int size){

    vector<vector<float>> M(size, vector<float>(size, 0));
    float sum=0;

    for (int i=0; i<size; i++){
        sum=0.0;
        for(int j=0; j<size; j++){
            M[i][j] = MIN_VALUE + static_cast<float>(rand()) * static_cast<float>(MAX_VALUE - MIN_VALUE) / RAND_MAX;
            sum+=M[i][j];
        }

        //strongly diagonal dominant
        M[i][i]=((float)2*(sum-M[i][i]));
    }
    return M;
}

vector<float> RHSVectorGenerator(int size){
    vector<float> b(size);
    for(int i=0; i<size; i++)
        b[i] = MIN_VALUE + static_cast<float>(rand()) * static_cast<float>(MAX_VALUE - MIN_VALUE) / RAND_MAX;

    return b;
}

void printResult(vector<float> b){
    for(int i=0; i<b.size(); i++)
            cout<< b[i] << endl;
}

void printMatrix(vector<vector<float>> M){
    for (int i=0; i<M[1].size(); i++){
        for(int j=0; j<M[1].size(); j++)
            cout<<M[i][j]<<" ";
        cout<<endl;
    }
}

#endif // UTILITIES_H
