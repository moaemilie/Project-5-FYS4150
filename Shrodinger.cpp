
#include <armadillo>
#include <iostream>
#include <math.h>
#include "Shrodinger.hpp"


Shrodinger::Shrodinger(double M_in, double h_in, double deltat_in){
    r = arma::cx_double(0.0, 1.0*deltat_in/(2*h_in*h_in));
    h = h_in;
    M = M_in;
    deltat = deltat_in;
}

double Shrodinger::getK(double i, double j, double N) {
    double k = j*N + i;
    return k;
}


arma::cx_mat Shrodinger::tri_matrix(arma::vec vector, arma::cx_double r){

    double N = vector.size();
    arma::cx_mat A = arma::cx_mat(N,N).fill(0.);

        //Create tridiagonal matrix
        A(0,0)=vector(0);
        A(0,1)=r;
        A(0,3)=r;
        A(N-1, N-1)=vector(N-1);

        for(int i=0; i < N-1; i++){
            for(int j=0;j<N-1; j++){

                 if(i==j){
                    A(i,j)=vector(i);

                    A(i+1,j) = r;
                    A(i,j+1) = r; 


                    if(j+3<N){
                        A(i+3,j) = r;
                        A(i,j+3) = r;
                    }

                    if((j+1)%3 == 0 && j!=N){
                        A(i+1,j) = 0;
                        A(i,j+1) = 0;
                    }
                
                }
        
            }
        }

    return A;
}

void Shrodinger::CalcAB(double M, double h, double deltat, arma::mat V){
    double k_len = (M-2)*(M-2);
    a_k = arma::cx_vec(k_len).fill(0.);
    b_k = arma::cx_vec(k_len).fill(0.);
    double n = V.size();
    //std::cout<< n;

    arma::cx_double i = arma::cx_double(0.0, 1.0);

    for(int i=0; i < M-2; i++){
            for(int j=0; j < M-2; j++){
            double k = getK(i, j, M-2);
            a_k(k) = arma::cx_double(1.0, 0.0) + arma::cx_double(4.0, 0.0)*r+ ((arma::cx_double(0.0, 1.0)*deltat)/(2.))*V(i,j);
            b_k(k) = arma::cx_double(1.0, 0.0) - arma::cx_double(4.0, 0.0)*r - ((arma::cx_double(0.0, 1.0)*deltat)/(2.))*V(i,j);
            }
        }


}