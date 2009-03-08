#include <src/StringMethod.h>

#include <iostream>

using namespace std;
using namespace ProtoMol;

int main(int argc, char *argv[]) {

   if (argc < 5) {
      cout <<"Usage : progname <string start filename> <threshold for commonspace> <do endpoints?> <tmdconfig>"<<endl;
      return 1;
   }

   int do_ends = atoi(argv[3]);

   StringMethod myString;

   myString.initialize(argv[1],argv[0],(double)atof(argv[2]), do_ends, argv[4]);

   while (myString.StringStep()) {

        myString.GetPhiPsiCoordinates();
        double *p = myString.ppp;

        myString.myPhiPsiOutput->run(p);

        continue;

   } 

   unsigned int sz = myString.apps.size();

   for(unsigned int i=0;i<sz;i++) myString.apps[i]->finalize();

   return 0;
}
