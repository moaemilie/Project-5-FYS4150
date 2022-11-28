
#include <armadillo>
#include <iostream>
#include <math.h>
#include "Shrodinger.cpp"
#include "Shrodinger.hpp"


int main(){

    //arma::cx_double r = arma::cx_double(0.0, 1.0);
    
    //std::complex<double>r i(0.0,1.0);
    arma::cx_vec vector = arma::cx_vec("0.0 1. 2. 3. 4. 5. 6. 7. 8.");
    //arma::vec vector = arma::vec("0.1 0.2 0.3 0.4 0.8 1. 5.");
    arma::mat V = arma::mat(3, 3).fill(1.);
    Shrodinger model = Shrodinger(0.2, 0.1, V, 1., 1., 1., 1., 1., 1.);
    arma::cx_mat A = model.tri_matrix(vector, -model.r);

    assert(A(0,0)==arma::cx_double(0.0, 0.0));
    assert(A(1,1)==arma::cx_double(1.0, 0.0));
    assert(A(2,2)==arma::cx_double(2.0, 0.0));
    assert(A(3,3)==arma::cx_double(3.0, 0.0));
    assert(A(4,4)==arma::cx_double(4.0, 0.0));
    assert(A(5,5)==arma::cx_double(5.0, 0.0));
    assert(A(6,6)==arma::cx_double(6.0, 0.0));
    assert(A(7,7)==arma::cx_double(7.0, 0.0));
    assert(A(8,8)==arma::cx_double(8.0, 0.0));

    assert(A(1,0)==-model.r);
    assert(A(2,1)==-model.r);
    assert(A(3,2)==arma::cx_double(0.0, 0.0));
    assert(A(4,3)==-model.r);
    assert(A(5,4)==-model.r);
    assert(A(6,5)==arma::cx_double(0.0, 0.0));
    assert(A(7,6)==-model.r);
    assert(A(7,8)==-model.r);

    assert(A(0,1)==-model.r);
    assert(A(1,2)==-model.r);
    assert(A(2,3)==arma::cx_double(0.0, 0.0));
    assert(A(3,4)==-model.r);
    assert(A(4,5)==-model.r);
    assert(A(5,6)==arma::cx_double(0.0, 0.0));
    assert(A(6,7)==-model.r);
    assert(A(8,7)==-model.r);

    assert(A(4,0)==arma::cx_double(0.0, 0.0));
    assert(A(0,4)==arma::cx_double(0.0, 0.0));

    //arma::cx_vec a = model.CalcAB(true);
    //arma::cx_vec b = model.CalcAB(false);

    //model.CalcAB(true);
    
    std::cout<<"\n\n a_k:";
    std::cout<< model.a_k;
    std::cout<<"\n\n b_k:";
    std::cout<< model.b_k;  
    std::cout<<"\n\n A:";
    std::cout<< model.A;  
    std::cout<<"\n\n B:";
    std::cout<< model.B;  


 

    return 0;
}