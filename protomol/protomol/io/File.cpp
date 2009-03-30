#include <protomol/io/File.h>

#include <protomol/base/SystemUtilities.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/Zap.h>
#include <protomol/type/String.h>

using namespace std;
using namespace ProtoMol;

#include <sstream>
#include <iostream>

//____ File
File::File() : mode(ios::in | ios::out) {}
File::File(std::ios::openmode mode) : mode(mode) {}

File::File(std::ios::openmode mode, const std::string &filename) :
  mode(mode), filename(filename) {
  open();
}

File::~File() {}

bool File::open() {
  if (is_open()) close();
  file.clear();

  try {
    file.open(filename.c_str(), mode);

    // Work around a boost::iostreams bug for F@H core
    if ((mode & ios::app) || (mode & ios::ate))
      file.seekg(0, ios::end);
    else file.seekg(0);

  } catch (const ios::failure &e) {}

  return is_open();
}

bool File::open(const string &filename) {
  setFilename(filename);
  return open();
}

bool File::open(const string &filename, ios::openmode mode) {
  setFilename(filename);
  this->mode = mode;
  return open();
}

void File::close() {
  if (is_open()) file.close();
}

bool File::is_open() {
  return file.is_open();
}

bool File::isAccessible() {
  return ProtoMol::isAccessible(filename);
}

void File::read(char *c, streamsize count) {
#ifdef __SUNPRO_CC
  //____ Sun WorkShop CC does not properly read more than one char ...
  for (streamsize i = 0; i < count; ++i)
    file.get(c[i]);

#else
  file.read(c, count);
#endif
}

string File::getline() {
  string res;
  bool ok = !file.fail() && !file.eof();

  std::getline(file, res);

  if (ok && !file.bad() && file.fail() && file.eof())
    file.clear(file.rdstate() & (~ios::failbit));

  return res;
}

unsigned int File::getLineTokens(vector<string> &tokens) {
  tokens.clear();

  stringstream ss(getline());
  string str;
  while (ss >> str)
    if (!str.empty())
      tokens.push_back(str);

  return tokens.size();
}
