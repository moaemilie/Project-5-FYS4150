
#include <armadillo>
#include <iostream>
#include <math.h>
#include "Shrodinger.hpp"


Shrodinger::Shrodinger(double h_in, double deltat_in, double x_c_in, double y_c_in, double sig_x_in, double sig_y_in, double p_y_in, double p_x_in, double v_0_in, double slit_width_in, double part_width_in, double wall_width_in, double x_pos_in, int slits_in){
    r = arma::cx_double(0.0, 1.0*deltat_in/(2*h_in*h_in));
    h = h_in; // Step size
    M = 1./(h_in)+1.; // Size of matrix
    deltat = deltat_in; // Time step
    // System configurations
    x_c = x_c_in;
    y_c = y_c_in;
    sig_x = sig_x_in;
    sig_y = sig_y_in;
    p_y = p_y_in; // momentum in y direction
    p_x = p_x_in; // momentum in x direction
    x_pos = x_pos_in;
    v_0 = v_0_in; // potential
    slit_width = slit_width_in;
    part_width =  part_width_in;
    wall_width = wall_width_in;
    slits = slits_in; // Number of slits in system
    V = init_V(); // Define the potantial wall 
    a_k = CalcAB(true); // Calculate the a_k vector
    b_k = CalcAB(false); // Calculate the b_k vector
    A = tri_matrix(a_k, -r); // Calculate the A matrix
    B = tri_matrix(b_k, r); // Calculate the B matrix
    u = init_u(); // // Calculate the inital u


}
// Function that returns the k position of a vector given i,j position in matrix with size N
double Shrodinger::getK(double i, double j, double N) {
    double k = j*N + i;
    return k;
}

// Function that creates the initial trididagnoal matrix A and B
arma::cx_mat Shrodinger::tri_matrix(arma::cx_vec vector, arma::cx_double r){
    double N = (M-2)*(M-2);
    double kmax = (M-2)*(M-2) - 1;
    double blocksize = (M-2);

    arma::cx_mat A = arma::cx_mat(N,N).fill(0.); // Define matrix that wil be returned

        for(int j=0; j < M-2; j++){
            for(int i=0; i < M-2; i++){

                    int k = getK(i, j, M-2);

                    A(k,k)=vector(k);

                    if(k < kmax and i < blocksize){
                        A(k+1,k) = r;
                        A(k,k+1) = r; 
                    }


                    if(k + blocksize <= kmax){
                        A(k+(M-2),k) = r;
                        A(k,k+(M-2)) = r;
                    }
                
                }
        
            }
        

    return A;
}

//Function that creates the ak and bk vectors used in matrix A and B
arma::cx_vec Shrodinger::CalcAB(bool aorb){
    double k_len = (M-2)*(M-2);
    // Create two vectores that wil store a_k and b_k
    a_k = arma::cx_vec(k_len).fill(0.);
    b_k = arma::cx_vec(k_len).fill(0.);

    // Run trough matrix to get values
    for(int j=0; j < M-2; j++){
            for(int i=0; i < M-2; i++){
            double k = getK(i, j, M-2);
            a_k(k) = arma::cx_double(1.0, 0.0) + arma::cx_double(4.0, 0.0)*r + ((arma::cx_double(0.0, 1.0)*deltat)/(2.))*V(i,j);
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

// Function creates the initial u vector
arma::cx_vec Shrodinger::init_u(){

    double k_len = (M-2)*(M-2);
    // Define vector that wil store u values
    arma::cx_vec u_0 = arma::cx_vec(k_len).fill(0.);
    
    // Runn trough matrix to calculate values
    for(int j=0; j < M-2; j++){
            for(int i=0; i < M-2; i++){
                double k = getK(i, j, M-2);
                double x = i*1/(M-2);
                double y = j*1/(M-2);
                u_0(k) =  exp(-pow(x-x_c, 2)/(2*pow(sig_x, 2)) - pow(y-y_c, 2)/(2*pow(sig_y, 2)) + arma::cx_double(0.0, 1.0)*p_x*(x-x_c) + arma::cx_double(0.0, 1.0)*p_y*(y-y_c));

            }
        }
    u_0 = normalise(u_0);
    return u_0;
}

// Function that creates the inital V matrix
arma::mat Shrodinger::init_V(){
    V = arma::mat(M-2, M-2).fill(0.);
    double walls_y = (1-slit_width*slits+part_width*(slits-1))/(2);
    double walls_x = (x_pos/h)-(wall_width/(2.*h));

    // V matrix looks different depending on number of slits

    if(slits == 3){

        for(double i = walls_x; i <=walls_x+(wall_width/h); i++){

         for(double j = 0; j < M-2; j++){

            if(j<= 75 || (85 <= j && j <= 95) || (105 <= j && j <= 115) || j >= 125){

                V(i, j) = v_0;
                }
            } 
        }
    }  

    if(slits == 2){

        for(double i = walls_x; i <=walls_x+(wall_width/h); i++){

         for(double j = 0; j < M-2; j++){

            if(j<= 85 || (95 <= j && j <= 105) || j >= 115){

                V(i, j) = v_0;
                }
            } 
        }
    }  

    if(slits == 1){

        for(double i = walls_x; i <=walls_x+(wall_width/h); i++){

         for(double j = 0; j < M-2; j++){

            if( j<= 99 || j >= 101){

                V(i, j) = v_0;
                }
            } 
        }

    }
return V;
    
}

// Function that finds the u for the next timestep
void Shrodinger::find_u_next(){

    arma::cx_vec b = B*u;
    u = arma::spsolve(A, b); 
}

// Function that calcualtes the total probability 
arma::cx_double Shrodinger::find_p_val(arma::cx_vec u){
    return sum(conj(u)%u); 

}

// Function that calcualtes the probability distribution 
arma::cx_vec Shrodinger::find_p_vec(arma::cx_vec u){
    return conj(u)%u; 

}