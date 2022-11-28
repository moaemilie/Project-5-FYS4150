#include <armadillo>
#include <iostream>
#include <math.h>
#include "Shrodinger.cpp"
#include "Shrodinger.hpp"

int main(){


    double h = 0.005;
    double deltatt = 2.5*pow(10,-5);
    double T = 0.008;
    double x_c = 0.25;
    double sig_x = 0.05;
    double p_x = 200;
    double y_c = 0.5;
    double sig_y = 0.05;
    double p_y = 0;
    double v_0 = 0;
    double slit_width = 0.05;
    double part_width = 0.05;
    double wall_width = 0.02;
    double x_pos = 0.5;

    Shrodinger model = Shrodinger(h, deltatt, x_c, y_c, sig_x, sig_y, p_y, p_x, v_0, slit_width, part_width, wall_width, x_pos);
    arma::cx_mat result= arma::cx_mat(pow(model.M-2,2), T/deltatt);

    for(double t = 0; t <= T; t+=deltatt){

        model.find_u_next();
        result.col(t/deltatt) = model.u;

    }
    
    std::string my_output_file_name =  "model1.csv";
    result.save(my_output_file_name, arma::csv_ascii);
    }