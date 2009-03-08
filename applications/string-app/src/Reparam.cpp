#include <src/Reparam.h>
#include <cmath>

using namespace std;

Reparam::Reparam(int npts, int cvs): numpoints(npts), cvsize(cvs) {
   stringPos.clear();
   S.clear();

   fq = 0;
   
}

Reparam::~Reparam() {
   int n = (int)S.size();
   for(int i=0;i<n;i++) {
      if (S[i] != NULL) delete [] S[i];
   }
   n = (int)initialStringPos.size();
   for(int i=0;i<n;i++) {
      if (initialStringPos[i] != NULL) delete [] initialStringPos[i];
   }

   n = (int)stringPos.size();

   for(int i=0;i<n;i++) {
     if (stringPos[i] != NULL) delete [] stringPos[i];
   }

   myFile.close();
}

void Reparam::initialize(int nn) {

   string ofilename("reparam.out");
   myFile.open(ofilename.c_str(),fstream::out);

   total_size = nn;
   
}

// q() function in Trevors' reparameterization code
int Reparam::qval(int p) {

   int tmpQ = 2;
    
   double myl = lval(p);
   // I need some diagnostics here
   //What kind?
   while (!((BigLVal(tmpQ-1) < myl) && (myl <= BigLVal(tmpQ)))) tmpQ += 1;
   
   return tmpQ;
     
}    

// l() function in Trevors' reparameterization code
double Reparam::lval(int p) {
  int R = numpoints;
  double bl = BigLVal(R);
  return (p-1)*bl/(R-1);


}

// L() function in Trevors' reparameterization code
double Reparam::BigLVal(int p) {

   double result = 0;
   if (p==1) return result;
   for(int j=2;j<p+1;j++) {
   /*for(int j=2;j<p;j++) {*/
      result += CollectiveVariableNorm(j-1);
   }
   return result;

}

double Reparam::CollectiveVariableNorm(int q) {
   double nm = 0;
   double x;

   for(int j=0;j<cvsize;j++) {
	   x = stringPos[q][j]-stringPos[q-1][j];
       nm += x*x;
   }      
   return sqrt(nm);
}


void Reparam::SetS() {

   for(int ii=0;ii<numpoints;ii++) {
     double *v = new double[total_size];
     S.push_back(v);
   }
}



void Reparam::doRun(double *a, double *b) {

  if (S.empty()) {
    SetS();

    for(int i=0;i<cvsize;i++) S[0][i] = a[i];

    for(int i=0;i<cvsize;i++) S[numpoints-1][i] = b[i];

  }

  for(int i=1;i<(numpoints-1);i++) {
    for(int j=0;j<cvsize;j++) {
       S[i][j] += b[j]/(numpoints-1);
    }
  }
    

}

void Reparam::run() {

  //int npts = stringPos.size();
  if (S.empty()) { 
    SetS();

  }

#if 1
   for(int i=0;i<numpoints;i++) {
      cout<<"Position of point "<<i<<" in C-space "<<endl;
      for(int j=0;j<cvsize;j++) cout<<stringPos[i][j]<<" ";
      cout<<endl;
   }
 

#endif  


    double *a = S[0];
    for(int i=0;i<cvsize;i++) a[i]=stringPos[0][i];
    
    a = S[numpoints-1];
    for(int i=0;i<cvsize;i++) a[i]=stringPos[numpoints-1][i];


  for(int i=2;i<numpoints;i++) {
     int myQ = qval(i);
     double rr = CollectiveVariableNorm(myQ-1);
     for(int jj=0;jj<cvsize;jj++) {
       S[i-1][jj] = stringPos[myQ-1-1][jj] +(lval(i) - BigLVal(myQ-1))*(stringPos[myQ-1][jj]-stringPos[myQ-1-1][jj])/rr;
     }

  }

    
   double str_norm = 0; 
   for(int i=1;i<numpoints;i++) {
      str_norm += Dist(S[i],S[i-1]);
   }

   cout<<"STRING LENGTH "<<str_norm<<endl;

   cout <<"++++++++++++++++ After reparameterization +++++++++++++++++++"<<endl;

   for(int i=0;i<numpoints;i++) {
       for(int j=0;j<cvsize;j++) {
         cout<<S[i][j]<<" ";
       }
       cout<<endl<<"************************"<<endl;
   }
       
}

void Reparam::SetInitialStringPositions() {

  if (S.empty()) {
    cerr << "String positions not set"<<endl;
    exit(2);
  }


  for(int i=0;i<numpoints;i++) {
      double *d = new double[cvsize];
      for(int j=0;j<cvsize;j++) {
          d[j] = S[i][j];
      }
      initialStringPos.push_back(d);
  }

}

double Reparam::Dist(double *x, double *y) {

   double nm = 0;

   for(int i=0;i<cvsize;i++) {
      double myd = (y[i]-x[i]);
      nm += myd*myd;
   }
   return sqrt(nm);
}

double Reparam::StringNorm() {

   if ((S.size() != initialStringPos.size()) && ((int)S.size() != numpoints)) {
      cerr << "Size mismatch in data structures holding initial and current string positions";
      exit(2);
   }

   //calculate the normed distance between the initial and final positions
   double s_norm = 0;

   for(int i=0;i<numpoints;i++) {
      s_norm += Dist(initialStringPos[i],S[i]);
   }

   return s_norm;  

}

void Reparam::initialize_stringPos_DataStructures(int nmp, int _3N) {

      for(int i=0;i<nmp;i++) {
       double *dv = new double[_3N];
       stringPos.push_back(dv);
   }


}

void Reparam::Straight_line_approximation() {

    if (S.empty()) {
    SetS();

  }

    double *a = S[0];
    for(int i=0;i<cvsize;i++) a[i]=stringPos[0][i];

    a = S[numpoints-1];
    for(int i=0;i<cvsize;i++) a[i]=stringPos[numpoints-1][i];

    for(int i=1;i<(numpoints-1);i++) {
       //double *ss = S[i];
       for(int j=0;j<cvsize;j++) {
           S[i][j] = S[i-1][j] + (S[numpoints-1][j]-S[0][j])/(numpoints-1);
       }
    }

}
