#include <protomol/base/SystemUtilities.h>

#include <protomol/base/StringUtilities.h>
#include <protomol/base/Exception.h>

#ifdef _WIN32

//____ Define the missing symbols from <unistd.h> for M$ ....
#include <direct.h>
#define CHDIR _chdir
#define PATHSEP '\\'
#define PATHSEPSTR "\\"
#define access _access
#include <fcntl.h>
#include <io.h>
#define F_OK 0
#define W_OK 2
#define R_OK 4

#else

#include <unistd.h>
#define CHDIR chdir
#define PATHSEP '/'
#define PATHSEPSTR "/"
#include <pwd.h>
#endif

#include <sys/stat.h>

#if defined (_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h> /* for GetFullPathName */

#else
#include <libgen.h>
#endif

using namespace std;

namespace ProtoMol {
//____ changeDirectory
  bool changeDirectory(const string &fileName) {
    char *confFile = (char *)fileName.c_str();
    char *currentdir = confFile;
    char *tmp = NULL;

#ifdef WIN32
    // Replace all '/' by '\'
    for (tmp = confFile; *tmp; ++tmp)
      if (*tmp == '/')
        *tmp = '\\';
#endif


    for (tmp = confFile; *tmp; ++tmp) ;

    // find final null
    for (; tmp != confFile && *tmp != PATHSEP; --tmp) ;

    // find last '/'
    if (tmp != confFile) {
      *tmp = 0;
      confFile = tmp + 1;
      if (CHDIR(currentdir))
        return false;
    } else if (*tmp == PATHSEP)  // config file in / is odd, but it might happen
      if (CHDIR(PATHSEPSTR))
        return false;

    return true;
  }

//____ isAccessible
  bool isAccessible(const string &fileName) {
    return ::access(fileName.c_str(), F_OK) == 0;
  }

//____ protomolAbort

  static void (*myAbortFunction)() = NULL;

  void protomolAbort() {
    if (myAbortFunction) (*myAbortFunction)();

    THROW("EXIT");
  }

//____ setProtomolAbort
  void setProtomolAbort(void (*abortFunction)()) {
    myAbortFunction = abortFunction;
  }

//____ protomolExit

  static void (*myExitFunction)() = NULL;

  void protomolExit() {
    if (myExitFunction) (*myExitFunction)();

    THROW("EXIT");
  }

//____ setProtomolExit
  void setProtomolExit(void (*exitFunction)()) {
    myExitFunction = exitFunction;
  }

//____ protomolStartSerial
  static void (*myStartSerial)(bool) = NULL;

  void protomolStartSerial(bool exludeMaster) {
    if (myStartSerial != NULL)
      (*myStartSerial)(exludeMaster);
  }

//____ setProtomolExit
  void setProtomolStartSerial(void (*startSerialFunction)(bool)) {
    myStartSerial = startSerialFunction;
  }

//____ protomolEndSerial
  static void (*myEndSerial)(bool) = NULL;

  void protomolEndSerial(bool exludeMaster) {
    if (myEndSerial != NULL)
      (*myEndSerial)(exludeMaster);
  }

//____ setProtomolExit
  void setProtomolEndSerial(void (*endSerialFunction)(bool)) {
    myEndSerial = endSerialFunction;
  }

//____ ISLITTLEENDIAN
  struct Endian {
    // Helper class to make sure that we get endianess correct ... M$
    static bool isLittleEndian() {
      unsigned int tmp = 1;
      return 0 != *(reinterpret_cast<const char *> ( &tmp));
    }
  };
  const bool ISLITTLEENDIAN = Endian::isLittleEndian();

//____ getUserName
  string getUserName() {
#ifdef _WIN32
    return "Win32";
#else
    if (getpwuid(getuid()) != NULL)
      return string(getpwuid(getuid())->pw_name);
    else
      return toString(getuid());
#endif
  }

  string getCanonicalPath(const string &path) {
    char buf[4096];

#ifdef _WIN32
    char *finalpart;
    DWORD len = GetFullPathName(path.c_str(), 4096, buf, &finalpart);
    if (len == 0 || len > 4095)
      THROW(string("GetFullPathName '") + path + "' failed.");

    return buf;

#else
    char tmp[path.length() + 3];

    // The file might not exist yet but its directory must.
    strcpy(tmp, path.c_str());
    char *dir = dirname(tmp);

    if (!realpath(dir, buf))
      THROW(string("realpath '") + path + "' failed.");

    strcpy(tmp, path.c_str());

    return string(buf) + "/" + basename(tmp);
#endif
  }
}
