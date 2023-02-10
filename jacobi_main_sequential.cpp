#include <iostream>
#include <stdlib.h>
#include <vector>

#include "sequentialJacobi.h"

using namespace std;

int main(int argc, char * argv[]){

    int n_iterations = 0;
    int dim_matrix = 0;
    int show_result = 0;

    //check if exist the first argument to set number of iterations or help guide
    if(argv[1] == NULL){
        n_iterations = 500;
        dim_matrix = 1000;
        show_result = 0;
    }
    else{
        if(argv[1] == "help" || argv[1][0] == 'H' || argv[1][0] == 'h'){
            cout<<"--- Help ---"<<endl;
            cout<<"./jacobi_seq n_iterations dim_matrix show_result"<<endl;
            cout<<"Parameters:"<<endl;
            cout<<"n_iterations: set number of iterations (DEFAULT: 500)"<<endl;
            cout<<"dim_matrix: set dimension of matrix (nxn) (DEFAULT: 1000)"<<endl;
            cout<<"show_result: set an integer flag to visualize the algorithm result (DEFAULT: 0)"<<endl;
            return 0;
        }
        else
            n_iterations = (atoi(argv[1]) < 1) ? 500 : atoi(argv[1]);
    }

    //check if exist the second argument to set dimension of matrix
    if(argv[2] == NULL){
        dim_matrix = 1000;
        show_result = 0;
    }
    else
        dim_matrix = (atoi(argv[2]) <= 1) ? 1000 : atoi(argv[2]);
    
    //check if exist the third argument to set flag to show the result of algorithm
    if(argv[3] == NULL)
        show_result = 0;
    else
        show_result = (atoi(argv[3]) != 0 && atoi(argv[3]) != 1) ? 0 : atoi(argv[3]);

    //generate matrix and vector random
    vector<vector<float>> A(dim_matrix, vector<float>(dim_matrix,0));
    vector<float> b(dim_matrix);
    A=matrixGenerator(dim_matrix);
    b=RHSVectorGenerator(dim_matrix);

    long time_seq;		//variable for sequence time
    int nIter = 0;		//variable to store number of iterations done

    cout<<"Sequential execution"<<endl;
    vector<float> resSeq=seqJacobi(n_iterations, dim_matrix, A, b, &time_seq, &nIter);

    cout<<"Time seq: "<<time_seq<<endl;

    if (show_result == 1){
	cout<<endl<<"Results: "<<endl;
	printResult(resSeq);
	cout<<"Computed with "<< nIter <<" iterations."<<endl;
    }

    return 0;
}
