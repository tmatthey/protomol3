/* -*- c++ -*- */
#ifndef REPARAM_H
#define REPARAM_H

#include <vector>
#include <iostream>
#include <fstream>

class Reparam {

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    public:
       Reparam(int n, int cv);

       ~Reparam();
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods in Reparam
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    void initialize(int nn);

    //void GetUpdatedStringPositions(nmInts, int cvs);

    //reparameterization class should only do the reparameterization calculation on the data
    // provided to itself. It should not have the responsibility of obtaining data.
    void run();
    double BigLVal(int p);
    double CollectiveVariableNorm(int q);
    double lval(int p);
    int qval(int p);

    void doRun(double *a, double *b);

    // Set initialStringPos so that we can calculate the distance of current string from 
    // starting string
    void SetInitialStringPositions();
    double StringNorm();

    void initialize_stringPos_DataStructures(int nmp, int _3N);

    void Straight_line_approximation();

   private:
     void SetS();
     double Dist(double *x, double *y);
    

    public:
    
       std::vector<double *> stringPos; //string positions from simulation
       std::vector<double *> S; //reparameterized string positions

       std::vector<double *> initialStringPos; /* initial string positions */

       int numpoints, cvsize; //size of collective variable set

       std::ofstream myFile;

       int fq;

       int total_size;

};

#endif /*REPARAM_H*/
