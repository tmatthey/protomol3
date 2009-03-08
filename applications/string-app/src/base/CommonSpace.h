/* -*- c++ -*- */
#ifndef COMMONSPACE_H
#define COMMONSPACE_H

#include <iostream>
#include <vector>

  class CommonSpace {

   
     //Constructors and destructors

    public:
       CommonSpace();

       ~CommonSpace();

    public:
       void initialize(int np, int _3n, int ss, double dp_t, double *v0); /* initializes data structures */
       void run(std::vector<double *>& eigenvectors) ; /* estimates the common space */

       

    private:
       void AddSpace(double *v);
       void vec_sub_mult(double dp, int vec_index); 
       void vec_normalize();

       double dotProduct(double *Vn, int i, int j);
       double array_norm();


    //data structures

    public:
      double *commonBasis; /* common basis vectors */
      int numpoints; /* no of intermediate string points */
      int basis_size; /* size of the common basis set */
      int subspace_size; /* no of low frequency eigenvectors at each intermediate point.*/
      int _3N ; /*size of each eigenvector */
     
    private:
      std::vector<double> dp; /* dot products */
      double *temp_vec;
      double dp_threshold;



  };

#endif /* COMMONSPACE_H */
