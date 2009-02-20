def try_dirs(paths):
    for parts in paths:
        path = ''
        for part in parts:
            path = os.path.join(path, part)
        if os.path.isdir(path): return path

    return None

def printenv(key):
    str = key + '='
    if os.environ.has_key(key): str += "'" + os.environ[key] + "'"
    print str

def boost_error(msg):
    print msg
    printenv('BOOST_HOME')
    printenv('BOOST_INCLUDE_PATH')
    printenv('BOOST_LIB_PATH')
    printenv('BOOST_LIB_SUFFIX')
    printenv('BOOST_VERSION')
    Exit(1)

def boost_check_version(context, version):
    context.Message("Checking for boost version %s..." % version)

    # Boost versions are in format major.minor.subminor
    v_arr = version.split(".")
    version_n = 0
    if len(v_arr) > 0: version_n += int(v_arr[0]) * 100000
    if len(v_arr) > 1: version_n += int(v_arr[1]) * 100
    if len(v_arr) > 2: version_n += int(v_arr[2])

    ret = context.TryRun("""#include <boost/version.hpp>
     int main() {return BOOST_VERSION == %d ? 0 : 1;}
     \n""" % version_n, '.cpp')[0]

    context.Result(ret)
    return ret


def boost_configure(conf, hdrs, libs, version = '1.37', lib_suffix = ''):
    env = conf.env

    boost_inc = None
    boost_lib = None
    boost_ver = version
    boost_lib_suffix = lib_suffix

    if os.environ.has_key('BOOST_VERSION'):
        boost_ver = os.environ['BOOST_VERSION']

    if os.environ.has_key('BOOST_HOME'):
        boost_home = os.environ['BOOST_HOME']

        path = try_dirs([[boost_home, 'include', 'boost'],
                         [boost_home, 'include',
                          'boost-' + boost_ver.replace('.', '_'), 'boost'],
                         [boost_home, 'boost']])
        if path:
            boost_inc = os.path.dirname(path)
        else:
            print "WARNING: No boost include path found in BOOST_HOME"

        path = os.path.join(boost_home, 'lib')
        if os.path.isdir(path):
            boost_lib = path
        else:
            print "WARNING: No boost lib path found in BOOST_HOME"


    if os.environ.has_key('BOOST_INCLUDE_PATH'):
        boost_inc = os.environ['BOOST_INCLUDE_PATH']

    if os.environ.has_key('BOOST_LIB_PATH'):
        boost_lib = os.environ['BOOST_LIB_PATH']

    if os.environ.has_key('BOOST_LIB_SUFFIX'):
        boost_lib_suffix = os.environ['BOOST_LIB_SUFFIX'].replace('.', '_')


    if boost_inc != None:
        env.Append(CPPPATH = [boost_inc])

    if boost_lib != None:
        env.Append(LIBPATH = [boost_lib])
        env.Append(RPATH = [boost_lib])


    # Check version
    if not conf.BoostVersion(boost_ver):
        boost_error("Wrong version")

    for name in hdrs:
        header = os.path.join('boost', name + '.hpp')
        if not conf.CheckCXXHeader(header):
            boost_error("boost header '" + name + ".hpp' not found.")

    for name in libs:
        libname = 'boost_' + name + boost_lib_suffix

        if not conf.CheckLib(libname):
            boost_error(libname + ' library not found.')
