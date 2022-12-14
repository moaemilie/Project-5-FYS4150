#include <armadillo>
#include <iostream>
#include <math.h>
#include "Shrodinger.cpp"
#include "Shrodinger.hpp"

#define xxx std::cout << "Im here:" << __FILE__ <<":" << __LINE__ << std::endl;

int main(){

    // Define the input values
    double h_in = 0.005;//0.05;
    double deltatt_in = 2.5*pow(10,-5);
    double T_in = 0.004;
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
    int slits_in = 1;

    // Create instance of model
    Shrodinger model = Shrodinger(h_in, deltatt_in, x_c_in, y_c_in, sig_x_in, sig_y_in, p_y_in, p_x_in, v_0_in, slit_width_in, part_width_in, wall_width_in, x_pos_in, slits_in);
 
    // Create a cube that wil contain the results
    arma::cx_cube result = arma::cx_cube((model.M), (model.M), T_in/deltatt_in).fill(arma::cx_double(0.0, 0.0));

    // Calculate values for t = 0
    arma::cx_vec p = conj(model.u)%model.u;
    arma::cx_mat matrix_p = arma::cx_mat((model.M-2), (model.M-2));
    
    // Make vector to matrix
    for(double coloumn = 0; coloumn < (model.M-2); coloumn++){
        for(double i = 0; i < (model.M-2); i++){
            matrix_p(i, coloumn) = p((model.M-2)*coloumn + i); 
            } 
        }

    // Create padded matrix
    arma::cx_mat padded_matrix_p = arma::cx_mat(model.M, model.M).fill(arma::cx_double(0.0, 0.0));
    // Insert matrix in padded matrix 
    padded_matrix_p.submat(1 , 1, model.M-2, model.M-2) = matrix_p;
    // Add matrix to cube
    result.slice(0) = padded_matrix_p;

    // Run trough every timestep and repeat 
    for(double t = deltatt_in; t <= T_in; t+=deltatt_in){
        model.find_u_next(); // Calculate the next u

        // Make vector to matrix 
        arma::cx_vec p = conj(model.u)%model.u;
        arma::cx_mat matrix_p = arma::cx_mat((model.M-2), (model.M-2));

       for(double coloumn = 0; coloumn < (model.M-2); coloumn++){
            for(double i = 0; i < (model.M-2); i++){

                matrix_p(i, coloumn) = p((model.M-2)*coloumn + i);   
                } 
            }
        
        // Add matrix to cube
        arma::cx_mat padded_matrix_new = arma::cx_mat(model.M, model.M).fill(arma::cx_double(0.0, 0.0));
        // Insert matrix in padded matrix 
        padded_matrix_new.submat(1 , 1, model.M-2, model.M-2) = matrix_p;
        // Add matrix to cube
        result.slice(t/deltatt_in) = padded_matrix_new;
    
        // Print status
        std::cout<< "\n\n" ;
        std::cout<< t/deltatt_in;
        std::cout<< " ";
        std::cout<< "out of";
        std::cout<< " ";
        std::cout<< T_in/deltatt_in; 
    }

    // Save result in .bin file
    result.save("Task9_test.bin");  
   
    return 0;
    }