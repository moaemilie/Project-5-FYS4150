
#include <armadillo>
#include <iostream>
#include <math.h>
#include "Shrodinger.hpp"


Shrodinger::Shrodinger(double h_in, double deltat_in, double x_c_in, double y_c_in, double sig_x_in, double sig_y_in, double p_y_in, double p_x_in, double v_0_in, double slit_width_in, double part_width_in, double wall_width_in, double x_pos_in){
    r = arma::cx_double(0.0, 1.0*deltat_in/(2*h_in*h_in));
    h = h_in;
    M = 1/h_in;
    deltat = deltat_in;
    V = init_V();
    a_k = CalcAB(true);
    b_k = CalcAB(false);
    B = tri_matrix(b_k, r);
    A = tri_matrix(a_k, -r);
    x_c = x_c_in;
    y_c = y_c_in;
    sig_x = sig_x_in;
    sig_y = sig_y_in;
    p_y = p_y_in;
    p_x = p_x_in;
    x_pos = x_pos_in;
    v_0 = v_0_in;
    slit_width = slit_width_in;
    part_width =  part_width_in;
    wall_width = wall_width_in;
    u = init_u();


}

double Shrodinger::getK(double i, double j, double N) {
    double k = j*N + i;
    return k;
}


arma::cx_mat Shrodinger::tri_matrix(arma::cx_vec vector, arma::cx_double r){

    double N = vector.size();
    arma::cx_mat A = arma::cx_mat(N,N).fill(0.);

        //Create tridiagonal matrix
        A(0,0)=vector(0);
        A(0,1)=r;
        A(0,3)=r;
        A(N-1, N-1)=vector(N-1);

        for(int i=0; i < N-1; i++){
            for(int j=0;j<N-1; j++){

                 if(i==j){
                    A(i,j)=vector(i);

                    A(i+1,j) = r;
                    A(i,j+1) = r; 


                    if(j+3<N){
                        A(i+3,j) = r;
                        A(i,j+3) = r;
                    }

                    if((j+1)%3 == 0 && j!=N){
                        A(i+1,j) = 0;
                        A(i,j+1) = 0;
                    }
                
                }
        
            }
        }

    return A;
}

//arma::cx_vec Shrodinger::CalcAB(double M, double h, double deltat, arma::mat V, bool aorb){
arma::cx_vec Shrodinger::CalcAB(bool aorb){
    double k_len = (M-2)*(M-2);
    a_k = arma::cx_vec(k_len).fill(0.);
    b_k = arma::cx_vec(k_len).fill(0.);
    double n = V.size();
    //std::cout<< n;

    arma::cx_double i = arma::cx_double(0.0, 1.0);

    for(int i=0; i < M-2; i++){
            for(int j=0; j < M-2; j++){
            double k = getK(i, j, M-2);
            a_k(k) = arma::cx_double(1.0, 0.0) + arma::cx_double(4.0, 0.0)*r+ ((arma::cx_double(0.0, 1.0)*deltat)/(2.))*V(i,j);
            b_k(k) = arma::cx_double(1.0, 0.0) - arma::cx_double(4.0, 0.0)*r - ((arma::cx_double(0.0, 1.0)*deltat)/(2.))*V(i,j);
            }
        }
    if(aorb){

        return a_k;
        }
    else{

        return b_k;
    } 

}

void Shrodinger::find_u_next(){

    arma::cx_vec b = B*u;
    u = solve(A, b);

}

arma::cx_vec Shrodinger::init_u(){

    double k_len = (M-2)*(M-2);
    arma::cx_vec u_0 = arma::cx_vec(k_len).fill(0.);

    for(int i=0; i < M-2; i++){
            for(int j=0; j < M-2; j++){
                double k = getK(i, j, M-2);
                double x = i*1/(M-2);
                double y = j*1/(M-2);
                u_0(k) =  exp(-pow(x-x_c, 2)/(2*pow(sig_x, 2))-pow(y-y_c, 2)/(2*pow(sig_y, 2))+arma::cx_double(0.0, 1.0)*p_x*(x-x_c)+arma::cx_double(0.0, 1.0)*p_y*(y-y_c));

            }
        }
    u_0 = normalise(u_0);
    return u_0;
}

arma::mat Shrodinger::init_V(){
    V = arma::mat(M-2, M-2).fill(0.);
    double walls_y = (1-slit_width*2+part_width)/(2);
    double walls_x = x_pos-(wall_width/(2.*h));


    for(double i = walls_x; i <=walls_x+(wall_width/h); i++){

        for(double j = 1; j < M/h; j++){

            if(j<= 85 || (95 < j && j < 105) || j > 115){

                V(i, j) = v_0;
            }
        }
    }
return V;
    
}