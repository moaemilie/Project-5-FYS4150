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


    arma::cx_cube result = arma::cx_cube((1/(h_in)-2), (1/(h_in)-2), 1).fill(arma::cx_double(0.0, 0.0));//T_in/deltatt_in
    //arma::cx_mat result= arma::cx_mat(pow(model.M-2,2), 4);
    //arma::cx_vec result= arma::cx_vec(T_in/deltatt_in+1);
    //arma::cx_vec result_1 = arma::cx_vec(pow(model.M-2,2));
    //arma::cx_vec result_2 = arma::cx_vec(pow(model.M-2,2));
    //arma::cx_vec result_3 = arma::cx_vec(pow(model.M-2,2));

    //arma::cx_vec result_1 = model.find_p_vec(model.u);
    //result.col(1) = conj(model.u)%model.u;//model.find_p_vec(model.u)
        // Make vector to matrix again!
    arma::cx_vec p = conj(model.u)%model.u;
    arma::cx_mat matrix_p = arma::cx_mat((1/(h_in)-2), (1/(h_in)-2));
    
    for(double coloumn = 0; coloumn < (1/(h_in)-2); coloumn++){
        for(double i = 0; i < (1/(h_in)-2); i++){

            matrix_p(i, coloumn) = p((1/(h_in)-2)*coloumn + i);   
            } 
        }
    //std::cout<< matrix_p;
    result.slice(0) = matrix_p;

/* 
    double counter = 0;


    for(double t = deltatt_in; t <= T_in; t+=deltatt_in){
        model.find_u_next();

        // Make vector to matrix again!
        arma::cx_vec p = conj(model.u)%model.u;
        arma::cx_mat matrix_p = arma::cx_mat((1/h_in-2), (1/h_in-2));
        double coloumn = 0;

        while(coloumn < 1/h_in-2){
            for(double i = 0; i < (1/h_in-2); i++){
                matrix_p(i, coloumn) = p((h_in-2)*coloumn + i);       
                } 
            }
            coloumn += 1;
        
        result.slice(t/deltatt_in) = matrix_p;
        


        //result.insert()
        //result(counter,0) = counter;
        counter += 1;
        //result(t/deltatt_in) = model.find_p_val(model.u);

        if(counter == 0.001/deltatt_in){
            result.col(2) = conj(model.u)%model.u;
            //result.col(1) = model.find_p_vec(model.u);
            //arma::cx_vec result_2 = model.find_p_vec(model.u);
        }
        else if(counter == 0.002/deltatt_in){
            result.col(3) = conj(model.u)%model.u;
            //arma::cx_vec result_3 = model.find_p_vec(model.u);
            //result.col(2) = model.find_p_vec(model.u);
        }
 */
/*         arma::cx_vec result_2 = model.u;

        std::cout<< "\n\n" ;
        std::cout<< t/deltatt_in;
        std::cout<< " ";
        std::cout<< "out of";
        std::cout<< " ";
        std::cout<< T_in/deltatt_in; 


}*/
    result.save("matrix_trial.bin");

 /*    std::string my_output_file_name =  "model6.csv";
    result.save(my_output_file_name, arma::csv_ascii);  */

    /* std::string filename = "test_model5.txt";
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
    ofile.close();   */

    return 0;
    }