#include <src/base/StringConfigReader.h>

using namespace std;

StringConfigReader::StringConfigReader() {}

StringConfigReader::StringConfigReader (const string &filename): Reader(filename) {}

StringConfigReader::~StringConfigReader(){}

bool StringConfigReader::tryFormat() {
   if (!open()) return false;
   else return true;
}


bool StringConfigReader::read() {
   if (!open()) return false;
    configfiles.clear();
   do {
        string fname(getline());
       std::cout<<fname<<std::endl;
       configfiles.push_back(fname);
    }while (!(file.eof()));
   configfiles.pop_back();
   std::cout<<configfiles.size()<<endl;
   for(unsigned int i=0;i<configfiles.size();i++) cout<<configfiles[i]<<endl;
   return true;
}
