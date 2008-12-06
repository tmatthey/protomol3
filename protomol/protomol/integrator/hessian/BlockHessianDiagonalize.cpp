#include <protomol/integrator/hessian/BlockHessianDiagonalize.h>
#include <protomol/integrator/hessian/BlockHessian.h>

#include <protomol/base/Report.h>

#include <protomol/type/BlockMatrix.h>

#if defined (HAVE_LAPACK)
#include <protomol/integrator/hessian/LapackProtomol.h>
#else
#if defined (HAVE_SIMTK_LAPACK)
#include "SimTKlapack.h"
#endif
#endif

#include <iostream>
#include <stdio.h>
#include <fstream>

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//Generate symetric or upper triangular Hessians
#define SYMHESS 1
//

namespace ProtoMol {
  //__________________________________________________ BlockHessianDiagonalize

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //constructors
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BlockHessianDiagonalize::BlockHessianDiagonalize() {
    bHess = 0; eigVal = 0; eigIndx = 0;
    blocks_num_eigs = 0; rE = 0;

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  BlockHessianDiagonalize::~BlockHessianDiagonalize() {		
    if(eigVal != 0)  delete [] eigVal;
    if(eigIndx != 0) delete [] eigIndx;
    if(blocks_num_eigs != 0) delete [] blocks_num_eigs;
    if(rE != 0) delete [] rE;

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Initialize for Block Hessians
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void BlockHessianDiagonalize::initialize(BlockHessian * bHessIn, const int sz) {

#if defined (HAVE_LAPACK)
#else
#if defined (HAVE_SIMTK_LAPACK)
#else
    THROW("Block Hessian diagonalization requires Lapack libraries.");
#endif
#endif
    //assign pointer to BlockHessian object    
    bHess = bHessIn;  
    //assign arrays
    try{
        //assign resudue eigenvector/value array	
        blocks_num_eigs = new int[bHess->num_blocks];
        rE = new double[bHess->hess_eig_size * 3];
        //
        eigVal = new double[sz];
        eigIndx = new int[sz];
    }catch(bad_alloc&){
        report << error << "[BlockHessianDiagonalize::initialize] Block Eigenvector array allocation error." << endr;
    }
    //Assign Eigenvector Blocks
    blockEigVect.resize(bHess->num_blocks);
    for(int i=0;i<bHess->num_blocks;i++){
      unsigned int start = bHess->hess_eig_point[i]*3;
      unsigned int rows = bHess->blocks_max[i]*3;
      blockEigVect[i].initialize(start,start,rows,rows);  //initialize block
    }
    //timers/counters for diagnostics
    rediagTime.reset();
    hessianTime.reset();

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Initialize for Full Hessians
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void BlockHessianDiagonalize::initialize(const int sz) {

#if defined (HAVE_LAPACK)
#else
#if defined (HAVE_SIMTK_LAPACK)
#else
    THROW("Block Hessian diagonalization requires Lapack libraries.");
#endif
#endif
    //assign pointer to BlockHessian object    
    bHess = 0;  
    //assign arrays
    try{
        //
        eigVal = new double[sz];
        eigIndx = new int[sz];
    }catch(bad_alloc&){
        report << error << "[BlockHessianDiagonalize::initialize] Eigenvector array allocation error." << endr;
    }
    //timers/counters for diagnostics
    rediagTime.reset();
    hessianTime.reset();

  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Find course eigenvectors, put into array with pointer mhQu, 
  // _3N rows, _rfM columns
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Real BlockHessianDiagonalize::findEigenvectors(const Vector3DBlock *myPositions,
                       const GenericTopology *myTopo, 
                       double * mhQu, const int sz_row, const int sz_col, 
                       const Real blockCutoffDistance, const Real eigenValueThresh) {

    //find 'minimum' Hessians for blocks
    hessianTime.start();	//time Hessian
    bHess->evaluateResidues(myPositions, myTopo, true);
    hessianTime.stop();	//stop timer
    //
    outputDiagnostics(3); //#### output Matrices
    //Diagonalize residue Hessians       
    //find coarse block eigenvectors 
    Real max_eigenvalue = findCoarseBlockEigs(eigenValueThresh);
    //
    report << hint << "[BlockHessianDiagonalize::findEigenvectors] Average residue eigenvalues "<<residues_total_eigs/bHess->num_blocks<<
              ", for "<<bHess->num_blocks<<" blocks."<<endr;
    //**** Find 'inner' Hessian H *********************************************
    //Find local block Hessian, and inter block Hessian within distance 'blockCutoffDistance'
    hessianTime.start();	//time Hessian
    bHess->clearBlocks();
    bHess->evaluateBlocks(blockCutoffDistance, myPositions, myTopo);
    hessianTime.stop();	//stop timer
    //put Q^T H Q into 'innerDiag' matrix
    innerHessian();
    //If cuttoff small do full calculation
    if(bHess->fullElectrostatics){
      //Q^T H Q for full electrostatics
      fullElectrostaticBlocks();
    }
    outputDiagnostics(0); //#### output Matrix
    //**** diagonalize 'inner' Hessian ****************************************
    int numeFound;
    rediagTime.start();    
    innerEigVec.initialize(0,0,residues_total_eigs,residues_total_eigs); //set small output matrix
    int info = diagHessian(innerEigVec.arrayPointer(), 
                            eigVal, innerDiag.arrayPointer(), residues_total_eigs, numeFound);
    rediagTime.stop();
    if( info == 0 ){            
      for(int i=0;i<residues_total_eigs;i++) eigIndx[i] = i;
      absSort(innerEigVec.arrayPointer(), eigVal, eigIndx, residues_total_eigs);
    }else{
      report << error << "[BlockHessianDiagonalize::findEigenvectors] Diagnostic diagonalization failed."<<endr;
    }
    //outputDiagnostics(1); //#### output eigenvalues
    //
    //**** Multiply block eigenvectors with inner eigenvectors ****************
    innerEigVec.columnResize(sz_col);
    for(int ii=0;ii<bHess->num_blocks;ii++)
      blockEigVect[ii].productToArray(innerEigVec, mhQu, 0, 0, sz_row, sz_col); //Aaa^{T}This
    //**** Return max eigenvalue for minimizer line search ********************
    return max_eigenvalue;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Find 'inner' matrix, S=Q^T H Q, using block matrix arithmetic
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void BlockHessianDiagonalize::innerHessian(){

    //Initialize storage
    innerDiag.initialize(0,0,residues_total_eigs,residues_total_eigs); //set small output matrix
    innerDiag.clear();
    //Re-size eigs space for num_eig vectors per block
    unsigned int i_res_sum = 0;
    for(int i=0;i<bHess->num_blocks;i++){
      blockEigVect[i].blockMove(blockEigVect[i].RowStart, i_res_sum);
      blockEigVect[i].columnResize(blocks_num_eigs[i]);
      i_res_sum += blocks_num_eigs[i];
    }
    //Do blocks if Hessian distance is 0
    for(int ii=0;ii<bHess->num_blocks;ii++){  
      BlockMatrix tempM((blockEigVect[ii]).ColumnStart, (bHess->blocks[ii]).ColumnStart, (blockEigVect[ii]).Columns, (bHess->blocks[ii]).Columns);
      (blockEigVect[ii]).transposeProduct(bHess->blocks[ii], tempM); //Aaa^{T}This
      tempM.product(blockEigVect[ii], innerDiag); //Aaa^{T}HAaa
    }
    //Do blocks if Hessian distance is 1
    for(int ii=0;ii<bHess->num_blocks-1;ii++){  
      BlockMatrix tempM((blockEigVect[ii]).ColumnStart, (bHess->adj_blocks[ii]).ColumnStart, (blockEigVect[ii]).Columns, (bHess->adj_blocks[ii]).Columns);
      (blockEigVect[ii]).transposeProduct(bHess->adj_blocks[ii], tempM); //Aaa^{T}This
      tempM.product(blockEigVect[ii+1], innerDiag); //Aaa^{T}HAbb
#if SYMHESS
      //Dont need this except for symmetric Hessians
      for(int jj=tempM.RowStart;jj<tempM.Rows+tempM.RowStart;jj++)
        for(int kk=(blockEigVect[ii+1]).ColumnStart;kk<(blockEigVect[ii+1]).Columns+(blockEigVect[ii+1]).ColumnStart;kk++)
           innerDiag(kk,jj) = innerDiag(jj,kk);
#endif
    }
    //Do non-adjacent bond blocks, distance > 1
    for(int ii=0;ii<bHess->non_adj_bond_blocks.size();ii++){  
      int ar0 = bHess->non_adj_bond_index[ii*2]; int ar1 = bHess->non_adj_bond_index[ii*2+1];
      BlockMatrix tempM((blockEigVect[ar0]).ColumnStart, (bHess->non_adj_bond_blocks[ii]).ColumnStart, (blockEigVect[ar0]).Columns, (bHess->non_adj_bond_blocks[ii]).Columns);
      (blockEigVect[ar0]).transposeProduct(bHess->non_adj_bond_blocks[ii], tempM); //Aaa^{T}This
      tempM.product(blockEigVect[ar1], innerDiag); //Aaa^{T}HAbb
#if SYMHESS
      //Dont need this except for symmetric Hessians
      for(int jj=tempM.RowStart;jj<tempM.Rows+tempM.RowStart;jj++)
        for(int kk=(blockEigVect[ar1]).ColumnStart;kk<(blockEigVect[ar1]).Columns+(blockEigVect[ar1]).ColumnStart;kk++)
          innerDiag(kk,jj) = innerDiag(jj,kk);
#endif
    }
    //Do adjacent non-bond blocks
    for(int ii=0;ii<bHess->adj_nonbond_blocks.size();ii++){  
      int ar0 = bHess->adj_nonbond_index[ii*2]; int ar1 = bHess->adj_nonbond_index[ii*2+1];
      BlockMatrix tempM((blockEigVect[ar0]).ColumnStart, (bHess->adj_nonbond_blocks[ii]).ColumnStart, (blockEigVect[ar0]).Columns, (bHess->adj_nonbond_blocks[ii]).Columns);
      (blockEigVect[ar0]).transposeProduct(bHess->adj_nonbond_blocks[ii], tempM); //Aaa^{T}This
      tempM.sumProduct(blockEigVect[ar1], innerDiag); //Aaa^{T}HAbb
#if SYMHESS 
      //Dont need this except for symmetric Hessians
      for(int jj=tempM.RowStart;jj<tempM.Rows+tempM.RowStart;jj++)
        for(int kk=(blockEigVect[ar1]).ColumnStart;kk<(blockEigVect[ar1]).Columns+(blockEigVect[ar1]).ColumnStart;kk++)
          innerDiag(kk,jj) = innerDiag(jj,kk);
#endif
    }
    //
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Find 'inner' matrix, S'=Q^T H' Q, for full electrostatic Hessian H'
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void BlockHessianDiagonalize::fullElectrostaticBlocks(){
  
    //Full electrostatics
    for(int ii=0;ii<bHess->num_blocks;ii++){
      BlockMatrix tempM((blockEigVect[ii]).ColumnStart, (bHess->electroStatics).ColumnStart, (blockEigVect[ii]).Columns, bHess->electroStatics.Columns);
      (blockEigVect[ii]).transposeProduct(bHess->electroStatics, tempM); //Aaa^{T}This
#if SYMHESS
      for(int ll=0;ll<bHess->num_blocks;ll++){// should be from ii unless symmetric Hessians
#else
      for(int ll=ii;ll<bHess->num_blocks;ll++){
#endif
        tempM.sumProduct(blockEigVect[ll], innerDiag); //Aaa^{T}HAbb
      }
    }
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Find isolated block eigenvectors
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Real BlockHessianDiagonalize::findCoarseBlockEigs(const Real eigenValueThresh){

    //size and position eigvects
    for(int i=0;i<bHess->num_blocks;i++){
      unsigned int start = blockEigVect[i].RowStart;
      blockEigVect[i].blockMove(start, start);  //make square
      blockEigVect[i].columnResize(blockEigVect[i].Rows);
    }
    //    
    residues_total_eigs = 0;	//find total residue modes to use
    Real max_eigenvalue = 0;
    for(int ii=0;ii<bHess->num_blocks;ii++){
      int numFound;
      rediagTime.start();
      int infor = diagHessian(blockEigVect[ii].arrayPointer(), &rE[bHess->hess_eig_point[ii] * 3], 
                                bHess->blocks[ii].arrayPointer(), bHess->blocks[ii].Rows, numFound);
      rediagTime.stop();
      if(infor) report << error << "[BlockHessianDiagonalize::findCoarseBlockEigs] Residue "<<ii+1<<" diagonalization failed."<<endr;
      for(int i=0;i<bHess->blocks_max[ii] * 3;i++) eigIndx[i] = i;
      absSort(blockEigVect[ii].arrayPointer(), &rE[bHess->hess_eig_point[ii] * 3], eigIndx, bHess->blocks_max[ii] * 3);
      //find number of eigs required
      for(int jj=0;jj<bHess->blocks_max[ii] * 3;jj++){
        if(fabs(rE[bHess->hess_eig_point[ii] * 3 + jj]) >= eigenValueThresh){
          blocks_num_eigs[ii] = jj;
          residues_total_eigs += jj;
          break;
        }
      }
      //
      Real tempf = rE[(bHess->hess_eig_point[ii] + bHess->blocks_max[ii] - 1) * 3];
      if(max_eigenvalue < tempf) max_eigenvalue = tempf;
      //
      //####for(int jj=0;jj<bHess->blocks_max[ii] * 3;jj++)
      //####  report << hint << "Loop "<<ii<<", Eig "<<jj<<", value "<<rE[bHess->hess_eig_point[ii]*3+jj]<<", index "<<eigIndx[jj]<<" eigs "<<blocks_num_eigs[ii]<<" tot "<<residues_total_eigs<<" Thresh "<<eigenValueThresh<<endr;
    }
    return max_eigenvalue;
    //
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Diagnostic output
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void BlockHessianDiagonalize::outputDiagnostics(int typ){    
    ofstream myFile;	//for diagnostic output

    if(typ == 0){
      //Output 'inner' matrix for diagnostics
      myFile.open("eigM", ofstream::out);
      myFile.precision(10);
      for(int jj=0;jj<residues_total_eigs*residues_total_eigs;jj++)
        myFile << jj / residues_total_eigs + 1 << " " << jj % residues_total_eigs + 1
                  << " " << innerDiag.MyArray[jj] << endl;
      myFile.close();
    }    
    if(typ == 1){
      //output eigenvalues for diagnostics
      myFile.open("eigRed", ofstream::out);
      myFile.precision(10);
      for(int jj=0;jj<residues_total_eigs;jj++)
        myFile << jj << "  " << eigVal[jj] << endl;
      //close file
      myFile.close();
    }
    if(typ == 3){
      //Output block Hessians
      myFile.open("blockH", ofstream::out);
      myFile.precision(10);    
      for(int ii=0;ii<bHess->num_blocks;ii++){
        int b_max = bHess->blocks_max[ii]*3; 
        int start_r = bHess->blocks[ii].RowStart;
        int start_c = bHess->blocks[ii].ColumnStart;
        for(int jj=0;jj<b_max*b_max;jj++){
          myFile << start_r + (jj / b_max) + 1<< " " << start_c + (jj % b_max) + 1
                  << " " << bHess->blocks[ii].MyArray[jj] << endl;
        }
      }
      myFile.close();
    }
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Diagonalize Hessian
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int BlockHessianDiagonalize::diagHessian(double *eigVecO, double *eigValO, double *hsnhessM, int dim, int &numFound){
   double *wrkSp;
   int *isuppz, *iwork;

   wrkSp = new double[26*dim];
   isuppz = new int[2*dim];
   iwork = new int[10*dim];
   //Diagonalize
#if defined(HAVE_LAPACK) || defined(HAVE_SIMTK_LAPACK)
    char *jobz = "V"; char *range = "A"; char *uplo = "U"; /* LAPACK checks only first character N/V */
    int n = dim;             /* order of coefficient matrix a  */
    int lda = dim;           /* leading dimension of a array*/
    double vl = 1.0;
    double vu = 1.0; 
    int il = 1;
    int iu = 1;
    double abstol = 0;
    int ldz = dim; int lwork = 26*dim; /* dimension of work array*///int m; 
    int liwork = 10*dim;						/* dimension of int work array*/
    //Recomended abstol for max precision
    char *cmach = "safe min";
#endif
    int info = 0;				/* output 0=success */
	int m = 0;
    //call LAPACK 
    //	
#if defined( HAVE_LAPACK )
    abstol = dlamch_( cmach);	//find machine safe minimum  
    //
    dsyevr_( jobz, range, uplo, &n, hsnhessM, &lda, &vl, &vu, &il, &iu, &abstol, &m, eigValO, eigVecO, &ldz, isuppz, 
                wrkSp, &lwork, iwork, &liwork, &info);
#else
#if defined( HAVE_SIMTK_LAPACK )
    int len_cmach = 8;
    int len_jobz = 1; int len_range = 1; int len_uplo = 1;
    abstol = dlamch_( *cmach, len_cmach);	//find machine safe minimum  
    //
    dsyevr_( *jobz, *range, *uplo, n, hsnhessM, lda, &vl, &vu, &il, &iu, &abstol, m, eigValO, eigVecO, ldz, isuppz, 
                wrkSp, lwork, iwork, &liwork, info, len_jobz, len_range, len_uplo);
#endif
#endif    
	numFound = m;
    //delete arrays
    delete [] iwork;
    delete [] isuppz;
    delete [] wrkSp;
    //return status
    return info;
  }

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Sort vectors for absolute value
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  void BlockHessianDiagonalize::absSort(double *eigVec, double *eigVal, int *eigIndx, int dim){
    int i;

    //find minimum abs value
    double minEv = fabs(eigVal[0]);
    for(i=1;i<dim;i++){
        if(minEv < fabs(eigVal[i])) break;
        else minEv = fabs(eigVal[i]);
    }
    i--;
    //sort around min
    if(i>0){
        int j = 0;
        eigIndx[j++] = i;
        int negp = i-1;
        int posp = i+1;
        while(negp >= 0 && posp < dim){
            if(fabs(eigVal[negp]) < fabs(eigVal[posp]))
                eigIndx[j++] = negp--;
            else eigIndx[j++] = posp++;
        }
        while(negp >= 0) eigIndx[j++] = negp--;
    }
    //Sort actual eigenvector array
    double *tmpVect = new double[dim];
    double tmpElt, tmpEval;
    int ii, k;
    for(i=0;i<dim;i++){
        if( eigIndx[i] != (int)i && eigIndx[i] != -1){		//need to swap?
            for(int j=0;j<dim;j++){
                tmpVect[j] = eigVec[i*dim+j];
                eigVec[i*dim+j] = eigVec[eigIndx[i]*dim+j];
            }
			//
			tmpEval = eigVal[i];
			eigVal[i] = eigVal[eigIndx[i]];
			//
            eigIndx[i] = -1;								//flag swapped
            ii = i;
            do{
                for(k=0;k<dim && eigIndx[k]!=(int)ii;k++);	//find where tmpVect goes
                if(k==dim || k==ii) break;					//end chain where indeces are equal
                for(int j=0;j<dim;j++){				//put it there
                    tmpElt = tmpVect[j];
                    tmpVect[j] = eigVec[k*dim+j];
                    eigVec[k*dim+j] = tmpElt;
                }
				//
				tmpElt = tmpEval;
				tmpEval = eigVal[k];
				eigVal[k] = tmpElt;
				//
                eigIndx[k] = -1;							//flag swapped
                ii = k;
            }while(k<dim);
        }
    }
    delete [] tmpVect;
  }

}
