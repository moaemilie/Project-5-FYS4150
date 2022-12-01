#include <armadillo>
#include <iostream>
#include <math.h>
#include "Shrodinger.cpp"
#include "Shrodinger.hpp"


int main(){


    double h_in = 0.05;
    double deltatt_in = 0.0001;//2.5*pow(10,-5); //This is the coccrect one!
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
    
    std::cout<< model.V;


}