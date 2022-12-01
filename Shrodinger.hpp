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
        arma::mat V;
        arma::cx_mat B;
        arma::cx_mat A;
        arma::cx_vec u;
        double x_c;
        double y_c;
        double sig_x; 
        double sig_y; 
        double p_y; 
        double p_x;
        double x_pos;
        double v_0;
        double slit_width;
        double part_width;
        double wall_width;



        Shrodinger(){};
        // Constructor
        Shrodinger(double h, double deltat, double x_c, double y_c, double sig_x, double sig_y, double p_y, double p_x, double v_0, double slit_width, double part_width, double wall_width, double x_pos);

        double getK(double i, double j, double N);

        arma::cx_mat tri_matrix(arma::cx_vec vector, arma::cx_double r);

        //arma::cx_vec CalcAB(double M, double h, double deltat, arma::mat V, bool aorb);
        arma::cx_vec CalcAB(bool aorb);

        arma::cx_vec init_u();

        arma::mat init_V();

        void find_u_next();

        arma::cx_double find_p_val(arma::cx_vec u);

        arma::cx_vec find_p_vec(arma::cx_vec u);
        

};

#endif