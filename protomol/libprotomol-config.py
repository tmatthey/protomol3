import os
from SCons.Script import *
import config
from platform import machine, architecture

deps = ['libfah']


def check_envvar(name, die = False):
    if os.environ.has_key(name):
        return os.environ[name]
    else:
        if die:
            print "Environment Variable Missing: " + name
            sys.exit(1)
        else: return None
        

def check_header(conf, name, die = False):
    if not conf.CheckCXXHeader(name):
        if die:
            print "Header Missing: " + name
            sys.exit(1)
        else: return False
    else: return True


def check_library(conf, name, die = False):
    if not conf.CheckLib(name):
        if die:
            print "Library Missing: " + name
            sys.exit(1)
        else: return False
    else: return True


def add_vars(vars):
    vars.AddVariables(
        BoolVariable('fah', 'Set to 1 to build library for Folding@Home', 0),
        BoolVariable('qrdiag', 'Set to 1 if QR diagonalization', 0),
        BoolVariable('gui', 'Set to 1 if using the GUI', 0),
        BoolVariable('lapack', 'Use LAPACK', 0),
        BoolVariable('simtk_lapack', 'Use SimTK LAPACK', 0),
        EnumVariable('openmm', 'Build with OpenMM', 'none',
                     allowed_values = ('none', 'reference', 'cuda')),
        BoolVariable('gromacs', 'Enable or disable gromacs support', 0),
        )


def CheckMKL(context):
    env = context.env

    context.Message('Checking for Intel MKL... ')

    save_CPPPATH = env['CPPPATH']
    save_LIBPATH = env['LIBPATH']

    if os.getenv('MKLROOT'):
        mkl_root = os.getenv('MKLROOT')
        env.Append(CPPPATH = [mkl_root + '/include'])
        env.Append(LIBPATH = [mkl_root + '/lib'])

    source = """
      #include "mkl.h"
      int main()
      {
        char ta = 'N', tb = 'N';
        int M = 1, N = 1, K = 1, lda = 1, ldb = 1, ldc = 1;
        double alpha = 1.0, beta = 1.0, *A = 0, *B = 0, *C = 0;
        dgemm(&ta, &tb, &M, &N, &K, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);
        return 0;
      }
    """
    source = source.strip()
    if architecture()[0] == '64bit': suffix = '_lp64'
    else: suffix = ''
    libs = ['mkl_intel' + suffix]

    # Thread model
    compiler = env.get('complier')
    if compiler == 'gnu': libs += ['mkl_gnu_thread', 'pthread']
    else: libs += ['mkl_intel_thread', 'iomp5']

    libs += ['mkl_core']

    if context.TryCompile(source, '.cpp'):
        save_LIBS = env['LIBS']

        env.Append(LIBS = libs)
        env.Append(LIBS = libs) # Twice to resolve all deps
        if context.TryLink(source, '.cpp'):
            env.Append(CPPDEFINES = ['HAVE_MKL_LAPACK'])
            context.Result(True)
            return True

        env.Replace(LIBS = save_LIBS)

    env.Replace(CPPPATH = save_CPPPATH)
    env.Replace(LIBPATH = save_LIBPATH)

    context.Result(False)
    return False


def configure_deps(conf):
    env = conf.env

    # libfah
    use_fah = int(env.get('fah', 0))
    if use_fah:
        config.configure('libfah', conf)

    # DIAG Options
    use_qr = int(env.get('qrdiag', 0))
    if use_qr == 1:
        env.Append(CCFLAGS = '-DHAVE_QRDIAG')

    # GUI Options
    use_gui = int(env.get('gui',0))
    if use_gui == 1:
        env.Append(CCFLAGS = '-DHAVE_GUI')

        if env['PLATFORM'] == 'win32':
            check_library(conf, 'wsock32', True)
        else:
            check_library(conf, 'pthread', True)

    # LAPACK
    have_lapack = False
    use_lapack = int(env.get('lapack', 1))
    if use_lapack:
        lapack_home = check_envvar('LAPACK_HOME')

        if lapack_home != None:
            env.Append(CPPPATH = [lapack_home])
            env.Append(LIBPATH = [lapack_home])

        if not hasattr(conf, 'CheckMKL'): conf.AddTest('CheckMKL', CheckMKL)
        have_lapack = conf.CheckMKL()

        if not have_lapack and (
            conf.CheckLib('lapack-3') or check_library(conf, 'lapack')):
            env.Append(CPPDEFINES = ['HAVE_LAPACK'])
            have_lapack = True

            # BLAS
            home = check_envvar('BLAS_HOME')
            if home != None: env.Append(LIBPATH = [home])
            conf.CheckLib('blas-3') or check_library(conf, 'blas')

            if env['PLATFORM'] in ['posix', 'darwin']:
                # G2C
                home = check_envvar('G2C_HOME')
                if home != None: env.Append(LIBPATH = [home])
                check_library(conf, 'g2c')

                # GFortran
                home = check_envvar('GFORTRAN_HOME')
                if home != None: env.Append(LIBPATH = [home])
                check_library(conf, 'gfortran')

    # SimTK LAPACK
    use_simtk = int(env.get('simtk_lapack', 0))
    if use_simtk == 1 or (use_lapack and not have_lapack):
        simtk_home = check_envvar('SIMTK_LAPACK_HOME')

        if simtk_home != None:
            env.Append(LIBPATH = [simtk_home + '/lib'])
            env.Append(CPPPATH = [simtk_home + '/include'])

        if check_library(conf, 'SimTKlapack') and \
                check_header(conf, 'SimTKlapack.h', True):
            env.Append(CPPDEFINES = ['HAVE_SIMTK_LAPACK'])
            have_lapack = True

    if not have_lapack:
        print "Missing lapack"
        sys.exit(1)

    # OpenMM Options
    openmm_type = env.get('openmm', 'none')
    if openmm_type != 'none':       
        # The following must bail if it is not found as openmm is not
        # installed to a place that the compiler will locate by default. The
        # same is also true for CUDA.
        openmm_home = check_envvar('OPENMM_HOME', True)

        env.Append(CPPPATH = [openmm_home + os.sep + 'include'])
        env.Append(LIBPATH = [openmm_home + os.sep + 'lib'    ])

        if openmm_type == 'reference':
            check_library(conf, 'OpenMM_d', True)
            env.Append(CPPDEFINES = ['HAVE_OPENMM'])

        if openmm_type == 'cuda':
            cuda_home = check_envvar('CUDA_HOME', True)

            env.Append(LIBPATH = [cuda_home + os.sep + 'lib'])

            check_library(conf, 'OpenMM_d', True)
            check_library(conf, 'OpenMMCuda_d', True)
            check_library(conf, 'cudart', True)

            env.Append(CPPDEFINES = ['HAVE_OPENMM'])

    # Gromacs
    gromacs = int(env.get('gromacs', 0))
    if gromacs:
        if os.environ.has_key('GROMACS_HOME'):
            ghome = os.environ['GROMACS_HOME']
            env.Append(LIBPATH = [ghome + '/lib'])
            env.Append(CPPPATH = [ghome + '/include'])

        check_library(conf, 'md', True)
        check_library(conf, 'gmx', True)

        for header in ['txtdump.h', 'names.h']:
            if not conf.CheckHeader('gromacs/' + header):
                print 'Need gromacs/' + header
                Exit(1)

        env.Append(CPPDEFINES = ['HAVE_GROMACS'])


def configure(conf):
    env = conf.env

    conf.AddTest('CheckMKL', CheckMKL)

    # Libprotomol
    if home:
        env.Append(CPPPATH = [home + '/src'])
        env.Append(LIBPATH = [home])

    lib = 'protomol'

    use_fah = int(env.get('fah', 0))
    if use_fah: lib = lib + '-fah'

    if not conf.CheckLib(lib):
        raise Exception, 'Need ' + lib + ' library';

    configure_deps(conf)
