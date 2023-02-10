#include <iostream>
#include <stdlib.h>
#include <vector>

#include "sequentialJacobi.h"
#include "fflowJacobi.h"

using namespace std;

int main(int argc, char * argv[]){

    int n_iterations = 0;
    int dim_matrix = 0;
    int n_threads = 0;
    int show_result = 0;

    //check if exist the first argument to set number of iterations
    if(argv[1] == NULL){
        n_iterations = 500;
        dim_matrix = 1000;
        n_threads = 2;
        show_result = 0;
    }
    else{
        if(argv[1] == "help" || argv[1][0] == 'H' || argv[1][0] == 'h'){
            cout<<"--- Help ---"<<endl;
            cout<<"./jacobi_ff n_iterations dim_matrix m_threads show_result"<<endl;
            cout<<"Parameters:"<<endl;
            cout<<"n_iterations: set number of iterations (DEFAULT: 500)"<<endl;
            cout<<"dim_matrix: set dimension of matrix (nxn) (DEFAULT: 1000)"<<endl;
            cout<<"n_threads: set number of thread (DEFAULT: 2)"<<endl;
            cout<<"show_result: set an integer flag to visualize the algorithm result (DEFAULT: 0)"<<endl;
            return 0;
        }
        else
            n_iterations = (atoi(argv[1]) < 1) ? 500 : atoi(argv[1]);
    }

    //check if exist the second argument to set dimension of matrix
    if(argv[2] == NULL){
        dim_matrix = 1000;
        n_threads = 2;
        show_result = 0;
    }
    else
        dim_matrix = (atoi(argv[2]) <= 1) ? 1000 : atoi(argv[2]);

    //check if exist the third argument to set number of threads
    if(argv[3] == NULL){
        n_threads = 2;
        show_result = 0;
    }
    else
        n_threads = (atoi(argv[3]) <= 1) ? 2 : atoi(argv[3]);
    
    //check if exist the last argument to set flag to show the result of algorithm
    if(argv[4] == NULL)
        show_result = 0;
    else
        show_result = (atoi(argv[4]) != 0 && atoi(argv[4]) != 1) ? 0 : atoi(argv[4]);

    //generate matrix and vector random
    vector<vector<float>> A(dim_matrix, vector<float>(dim_matrix,0));
    vector<float> b(dim_matrix);
    A=matrixGenerator(dim_matrix);
    b=RHSVectorGenerator(dim_matrix);

    long time_seq;                  //variable for sequence time
    long time_ffN;                  //variable for fastflow time (n_threads>1)
    long time_ff1;                  //variable for fastflow time (n_threads=1)
    int nIter = n_iterations;

    vector<float> resSeq=seqJacobi(n_iterations, dim_matrix, A, b, &time_seq, &nIter);

    cout<<"Parallel execution with FastFlow"<<endl;
    vector<float> resFF=fflowJacobi(n_iterations, dim_matrix, n_threads, A, b, &time_ffN);
    vector<float> resFF1=fflowJacobi(n_iterations, dim_matrix, 1, A, b, &time_ff1);

    float speedUpFF=speedup(time_seq, time_ffN);
    float scalFF=scalability(time_ff1, time_ffN);
    float effFF=efficiency(time_seq, time_ffN, n_threads);

    cout<<"Speedup: "<<speedUpFF<<endl;
    cout<<"Scalability: "<<scalFF<<endl;
    cout<<"Efficiency: "<<effFF<<endl;

    if (show_result == 1){
    	cout<<endl<<"Results: "<<endl;
    	printResult(resSeq);
    }
        
    return 0;
}
