#include <src/base/StringPhiPsiOutput.h>
#include <iostream>

using namespace std;

StringPhiPsiOutput::StringPhiPsiOutput() : MyOutput() {}

StringPhiPsiOutput::StringPhiPsiOutput(string fname, int n) : MyOutput(fname), numpoints(n) {}

StringPhiPsiOutput::~StringPhiPsiOutput() {}

void StringPhiPsiOutput::doRun() {

   if (!myfile.is_open()) {
      cerr << "Output file is not open";
      exit(1);
   }

   for(int i=0;i<numpoints;i++) {
     myfile << myData[i*2+0]<<" " ; //phi value
     myfile << myData[i*2+1]<<" "; //psi value
   }
   myfile << endl;

}
