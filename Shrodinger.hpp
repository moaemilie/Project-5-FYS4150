#ifndef __filename_hpp__
#define __filename_hpp__
#include <armadillo>
#include <iostream>

class Shrodinger{

    private:

    public:
        arma::cx_double r;
        double h; // Step size
        double M; // Matrix size
        double deltat; // time step
        int slits; // Number of slits
        arma::cx_vec a_k; // vector with a_k values
        arma::cx_vec b_k; // vector with b_k values
        arma::mat V; // Matrix that stores the potential 
        arma::sp_cx_mat B; // B matrix
        arma::sp_cx_mat A; // A matrix
        arma::cx_vec u; // Vector that stores the u solution for the system in this instance
        // Simulation configurations:
        double x_c; 
        double y_c;
        double sig_x; 
        double sig_y; 
        double p_y; // momentum in y direction
        double p_x; // momentum in x direction
        double x_pos;
        double v_0; // Potential
        double slit_width; 
        double part_width;
        double wall_width;



        Shrodinger(){};
        // Constructor
        Shrodinger(double h, double deltat, double x_c, double y_c, double sig_x, double sig_y, double p_y, double p_x, double v_0, double slit_width, double part_width, double wall_width, double x_pos, int slits_in);

        // Function that returns the k position of a vector given i,j position in matrix with size N
        double getK(double i, double j, double N);

        // Function that creates the initial trididagnoal matrix A and B
        arma::cx_mat tri_matrix(arma::cx_vec vector, arma::cx_double r);

        //Function that creates the ak and bk vectors used in matrix A and B
        arma::cx_vec CalcAB(bool aorb);

        // Function creates the initial u vector
        arma::cx_vec init_u();

        // Function that creates the inital V matrix
        arma::mat init_V();

        // Function that finds the u for the next timestep
        void find_u_next();

        // Function that calcualtes the total probability 
        arma::cx_double find_p_val(arma::cx_vec u);

        // Function that calcualtes the probability distribution 
        arma::cx_vec find_p_vec(arma::cx_vec u);
        

};

#endif