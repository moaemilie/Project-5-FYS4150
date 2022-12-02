#include <armadillo>
#include <iostream>
#include <math.h>
#include "Shrodinger.cpp"
#include "Shrodinger.hpp"

int main(){


    double h_in = 0.005;
    double deltatt_in = 2.5*pow(10,-5); //This is the coccrect one!
    double T_in = 0.002;
    double x_c_in = 0.25;
    double sig_x_in = 0.05;
    double p_x_in = 200;
    double y_c_in = 0.5;
    double sig_y_in = 0.20;
    double p_y_in = 0;
    double v_0_in = 1*pow(10,10);
    double slit_width_in = 0.05;
    double part_width_in = 0.05;
    double wall_width_in = 0.02;
    double x_pos_in = 0.5;

    Shrodinger model = Shrodinger(h_in, deltatt_in, x_c_in, y_c_in, sig_x_in, sig_y_in, p_y_in, p_x_in, v_0_in, slit_width_in, part_width_in, wall_width_in, x_pos_in);
 
    // Create a cube that wil hold the results
    arma::cx_cube result = arma::cx_cube((1/(h_in)-2), (1/(h_in)-2), T_in/deltatt_in).fill(arma::cx_double(0.0, 0.0));

    // Calculate values for t = 0
    arma::cx_vec p = conj(model.u)%model.u;
    arma::cx_mat matrix_p = arma::cx_mat((1/(h_in)-2), (1/(h_in)-2));
    
    // Make vector to matrix again and add to cube
    for(double coloumn = 0; coloumn < (1/(h_in)-2); coloumn++){
        for(double i = 0; i < (1/(h_in)-2); i++){

            matrix_p(i, coloumn) = p((1/(h_in)-2)*coloumn + i);   
            } 
        }
    result.slice(0) = matrix_p;


    for(double t = deltatt_in; t <= T_in; t+=deltatt_in){
        model.find_u_next();

        // Make vector to matrix again
        arma::cx_vec p = conj(model.u)%model.u;
        arma::cx_mat matrix_p = arma::cx_mat((1/h_in-2), (1/h_in-2));
        for(double coloumn = 0; coloumn < (1/(h_in)-2); coloumn++){
            for(double i = 0; i < (1/(h_in)-2); i++){

                matrix_p(i, coloumn) = p((1/(h_in)-2)*coloumn + i);   
                } 
            }
        
        // Add matrix to cube
        result.slice(t/deltatt_in) = matrix_p;
 

        std::cout<< "\n\n" ;
        std::cout<< t/deltatt_in;
        std::cout<< " ";
        std::cout<< "out of";
        std::cout<< " ";
        std::cout<< T_in/deltatt_in; 
    }
    
    result.save("Task7_trial2_t_0.004.bin");


    return 0;
    }