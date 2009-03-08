#include <src/base/CommonSpace.h>
#include <cmath>

//#define DEBUG_COMMONSPACE

using namespace std;


CommonSpace::CommonSpace() {
  commonBasis = NULL;
  temp_vec = NULL;
}

CommonSpace::~CommonSpace() {

  if (commonBasis != NULL) delete [] commonBasis;

  if (temp_vec != NULL) delete [] temp_vec;
}

void CommonSpace::initialize(int np, int _3n, int ss, double dp_t, double *v0) {

   numpoints = np;

   _3N = _3n;

   subspace_size = ss;

   basis_size = 0;

   //initialize common basis size to _3Nx_3N
   if (commonBasis == NULL) commonBasis = new double[_3N*_3N];

   if (temp_vec == NULL) temp_vec = new double[_3N];

   //copy the first subspace_size vectors into the basis...they will always be members
   //of basis set.
   //memcpy(commonBasis, v0, subspace_size*_3N*sizeof(double));
   for(int i=0;i<_3N*subspace_size;i++) commonBasis[i] = v0[i];
   basis_size = subspace_size;

   dp_threshold = dp_t;

#ifdef DEBUG_COMMONSPACE
   //for (int k=0;k<subspace_size*_3N;k++) cout<<v0[k]<<" "<<commonBasis[k]<<endl;
#endif

}


void CommonSpace::run(vector<double *>& eigenvectors) {

   basis_size = subspace_size;

   unsigned int sz = eigenvectors.size();

   cout<<" Eigenvector structures/points "<<sz<<endl;
   
   for (unsigned int j=1;j<sz;j++) {

#if 0 
      double **Vv = eigenvectors[j]; 
      cout<<"Eigenvectors for point "<<j<<endl;
      cout<<(*Vv)[0]<<" "<<(*Vv)[_3N]<<" "<<(*Vv)[_3N*2]<<endl;
#endif

     AddSpace(eigenvectors[j]);

   }

   cout<<"Basis size "<<basis_size<<", Subspace size "<<subspace_size<<endl;
#ifdef DEBUG_COMMONSPACE

   for(int i=0;i<basis_size;i++)
      for(int k=0;k<_3N;k++) 
         cout<<commonBasis[i*_3N+k]<<endl;
#endif

}  

void CommonSpace::AddSpace(double *Vn) {


  for (int i=0;i<subspace_size;i++) {
     
     dp.clear(); /* reset dot product vector */

     for(int j=0;j<basis_size;j++) {

       dp.push_back(dotProduct(Vn,i,j)); 

     }

     double dp_norm = array_norm();

     cout<<"i "<<i<<", dp_norm "<<dp_norm<<endl;

     if (dp_norm > dp_threshold) {
        //ith vector is already contained in the subspace
        continue;
     } else {

        cout<<" i "<<i<<endl;

        //memcpy(temp_vec,Vn+i*_3N,_3N*sizeof(double)); /* copying ith eigenvector in temp_vec */
        for(int kk=0;kk<_3N;kk++) {
           temp_vec[kk] = Vn[i*_3N+kk];
        } /* getting rid of memcpy */
        for(int j=0;j<basis_size;j++) {
           vec_sub_mult(dp[j],j);
        }

        /* normalize temp_vec */
        vec_normalize();

        //memcpy(commonBasis+basis_size*_3N,temp_vec,_3N*sizeof(double));
        for(int kk=0;kk<_3N;kk++) {
           commonBasis[basis_size*_3N+kk] = temp_vec[kk];
        }
        basis_size++;

     }

  }

} 

double CommonSpace::dotProduct(double *Vn, int i, int j) {

   double dp = 0;

   //int v0_index = (j-1)*_3N;
   int v0_index = j*_3N;

   //int vn_index = (i-1)*_3N;
   int vn_index = i*_3N;

   for(int k=0;k<_3N;k++) {
#ifdef DEBUG_COMMONSPACE
      cout<<"Indices "<<v0_index+k<<", "<<vn_index+k<<endl;
      cout<<"commonBasis "<<commonBasis[v0_index+k]<<" ";
      cout<<"Vn "<<(*Vn)[vn_index+k]<<endl;
#endif

      dp += commonBasis[v0_index+k]*Vn[vn_index+k];
   }

   return dp;

}

double CommonSpace::array_norm() {

   double dp_sum=0;

   for(unsigned int i=0;i<dp.size();i++) dp_sum += dp[i]*dp[i];

   return sqrt(dp_sum);

}

void CommonSpace::vec_sub_mult(double dp, int vec_index) {

   for(int i=0;i<_3N;i++) temp_vec[i] = temp_vec[i] - dp*(commonBasis[(vec_index-1)*_3N+i]);

}

void CommonSpace::vec_normalize() {

   double nm = 0;

   for (int i=0;i<_3N;i++) nm += temp_vec[i]*temp_vec[i];

   nm = sqrt(nm);

   for (int i=0;i<_3N;i++) temp_vec[i] /= nm;

}
