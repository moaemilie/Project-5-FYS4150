#include <armadillo>
#include <iostream>
#include <math.h>
#include "Shrodinger.cpp"
#include "Shrodinger.hpp"

int main(){


    double h_in = 0.005;
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
 /*    
    model.find_u_next();

    std::cout<< "\n\n V dim:";
    std::cout<< size(model.V);
    std::cout<< "\n\n A dim:";
    std::cout<< size(model.A);
    std::cout<< "\n\n B dim:";
    std::cout<< size(model.B);
    std::cout<< "\n\n b dim:";
    std::cout<< size(model.b_k);
    std::cout<< "\n\n u dim:";
    std::cout<< size(model.u); */



    //arma::cx_mat result= arma::cx_mat(pow(model.M-2,2), 3);
    //arma::cx_vec result= arma::cx_vec(T_in/deltatt_in+1);
    //arma::cx_vec result_1 = arma::cx_vec(pow(model.M-2,2));
    arma::cx_vec result_2 = arma::cx_vec(pow(model.M-2,2));
    arma::cx_vec result_3 = arma::cx_vec(pow(model.M-2,2));

    arma::cx_vec result_1 = model.find_p_vec(model.u);
    //result.col(0) = model.find_p_vec(model.u);
    
    for(double t = deltatt_in; t <= T_in; t+=deltatt_in){
        model.find_u_next();
        //result(t/deltatt_in) = model.find_p_val(model.u);

        if(t == 0.001){
            //result.col(1) = model.find_p_vec(model.u);
            arma::cx_vec result_2 = model.find_p_vec(model.u);
        }
        else if(t == 0.002){
            //result.col(2) = model.find_p_vec(model.u);
            arma::cx_vec result_3 = model.find_p_vec(model.u);
        }

        std::cout<< "\n\n" ;
        std::cout<< t/deltatt_in;
        std::cout<< " ";
        std::cout<< "out of";
        std::cout<< " ";
        std::cout<< T_in/deltatt_in;

    } 

    //std::string my_output_file_name =  "model5.csv";
    //result.save(my_output_file_name, arma::csv_ascii); 

    std::string filename = "test_model5.txt";
    std::ofstream ofile;
    ofile.open(filename);
    int width = 12;
    int prec  = 4;
    // Loop over steps
    for (int i = 0; i < result_3.size(); i++)
    {
    ofile << std::setw(width) << std::setprecision(prec) << std::scientific << result_1[i]
            << std::setw(width) << std::setprecision(prec) << std::scientific << result_2[i]
            << std::setw(width) << std::setprecision(prec) << std::scientific << result_3[i]
            << std::endl;
    }  
    ofile.close();  

    return 0;
    }