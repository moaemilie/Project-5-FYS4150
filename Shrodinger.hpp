#ifndef __filename_hpp__
#define __filename_hpp__
#include <armadillo>
#include <iostream>

class Shrodinger{

    private:

    public:
        arma::cx_double r;
        double h;
        double M;
        double deltat;
        arma::cx_vec a_k;
        arma::cx_vec b_k;


        Shrodinger(){};
        // Constructor
        Shrodinger(double M, double h, double deltat);

        double getK(double i, double j, double N);

        arma::cx_mat tri_matrix(arma::vec vector, arma::cx_double r);

        void CalcAB(double M, double h, double deltat, arma::mat V);

};

#endif