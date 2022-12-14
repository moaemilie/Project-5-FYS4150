#include <armadillo>
#include <iostream>
#include <math.h>
#include "Shrodinger.cpp"
#include "Shrodinger.hpp"

int main(){

    // Define initial values
    double h_in = 0.005;
    double deltatt_in = 2.5*pow(10,-5); 
    double T_in = 0.008;
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
    int slits_in = 2;

    // Create instance of model
    Shrodinger model = Shrodinger(h_in, deltatt_in, x_c_in, y_c_in, sig_x_in, sig_y_in, p_y_in, p_x_in, v_0_in, slit_width_in, part_width_in, wall_width_in, x_pos_in, slits_in);

    // Create vector to store result
    arma::cx_vec result = arma::cx_vec(T_in/deltatt_in+1);

    result(0) = model.find_p_val(model.u); // Fill with inital probability

    // Run for every timestep
    for(double t = deltatt_in; t <= T_in; t+=deltatt_in){
        
        model.find_u_next(); //Calculate next u
        result(t/deltatt_in) = model.find_p_val(model.u); // Add probability to result vector
 
        // print status
        std::cout<< "\n\n" ;
        std::cout<< t/deltatt_in;
        std::cout<< " ";
        std::cout<< "out of";
        std::cout<< " ";
        std::cout<< T_in/deltatt_in; 
}

    // Save result in a .csv
    std::string my_output_file_name =  "Task7.csv";
    result.save(my_output_file_name, arma::csv_ascii);  

    return 0;
    }