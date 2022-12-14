# Project-5-FYS4150


- Shrodinger.hpp 

  We define the model that contains our Shcrodinger equation solver

- Shrodinger.cpp 

  We implement the functions for the class Shcrodinger and also initializes the attributes of the class.

- Task7.cpp 

  We make an instance of the Shcrodinger class with initial values and a double slitt barier. It then proseeds to run the simulation for T_in/deltatt_in timesteps.
  For every timestep the sum of the probability is stored in an arma cx vector. This is the writen to a .csv file.
  To simulate with and without a double-slit barrier by changig the input v_0_in. 

  Its compiled by:
  g++ -c Task7.cpp -std=c++11
  g++ Task7.o -o Task7.exe -larmadillo

  and run by:
  ./Task7.exe


- Task8.cpp 

  We make an instance of the Shcrodinger class with initial values and a double slitt barier. It then proseeds to run the simulation for T_in/deltatt_in  timesteps.
  For every timestep the probability is stored in an arma cx matrix which again is stored in a arma cx cube. This is the writen to a .bin file.

  Its compiled by:
  g++ -c Task8.cpp -std=c++11
  g++ Task8.o -o Task8.exe -larmadillo

  and run by:
  ./Task8.exe

- Task9.cpp

  We make an instance of the Shcrodinger class with initial values and a adjustable number of slits. It then proseeds to run the simulation for T_in/deltatt_in timesteps.
  For every timestep the probability is stored in an arma cx matrix which again is stored in a arma cx cube. This is the writen to a .bin file.

  Its compiled by:
  g++ -c Task9.cpp -std=c++11
  g++ Task9.o -o Task9.exe -larmadillo

  and run by:
  ./Task9.exe


- TestTriDiag.cpp
  
  This file test the function tri_matrix() in the Schrodinger class. 
  This was only used for testing for one spesific configuration. Its therefore no longer possible to run the code afther we have generalized it.
