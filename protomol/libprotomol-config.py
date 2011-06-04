import os
from SCons.Script import *
import config
from platform import machine, architecture

deps = ['libfah', 'mkl', 'lapack']


def add_vars(vars):
    vars.AddVariables(
        BoolVariable('fah', 'Set to 1 to build library for Folding@Home', 0),
        BoolVariable('qrdiag', 'Set to 1 if QR diagonalization', 0),
        BoolVariable('gui', 'Set to 1 if using the GUI', 0),
        EnumVariable('lapack', 'Use LAPACK', 'any', allowed_values =
                     ('1', 'any', 'none', 'mkl', 'simtk', 'system')),
        EnumVariable('openmm', 'Build with OpenMM', 'none',
                     allowed_values = ('none', 'reference', 'cuda')),
        BoolVariable('gromacs', 'Enable or disable gromacs support', 0),
        )


def configure_deps(conf):
    env = conf.env

    # libfah
    if env.get('fah', 0): config.configure('libfah', conf)

    # DIAG Options
    if env.get('qrdiag', 0): env.AppendUnique(CPPDEFINES = ['HAVE_QRDIAG'])

    # GUI Options
    if env.get('gui',0):
        if env['PLATFORM'] == 'win32': config.require_lib(conf, 'wsock32')
        else: config.require_lib(conf, 'pthread')

        env.AppendUnique(CPPDEFINES = ['HAVE_GUI'])

    # LAPACK
    have_lapack = False
    lapack = env.get('lapack', 'any')
    if lapack == '1' or lapack is True: lapack = 'any'
    elif lapack == '0' or lapack is False: lapack = 'none'

    if lapack != 'none':
        # Intel MKL LAPACK
        if not have_lapack and lapack in ['any', 'mkl']:
            have_lapack = config.configure('mkl', conf)

        # System LAPACK
        if not have_lapack and lapack in ['any', 'system']:
            have_lapack = config.configure('lapack', conf)

        # SimTK LAPACK
        if not have_lapack and lapack in ['any', 'simtk']:
            config.check_home(conf, 'simtk_lapack')

            if (config.check_lib(conf, 'SimTKlapack') and
                config.check_cxx_header(conf, 'SimTKlapack.h')):

                env.AppendUnique(CPPDEFINES = ['HAVE_SIMTK_LAPACK'])
                have_lapack = True

        if not have_lapack: raise Exception, "Missing LAPACK"


    # OpenMM Options
    openmm_type = env.get('openmm', 'none')
    if openmm_type != 'none':
        # The following must bail if it is not found as openmm is not
        # installed to a place that the compiler will locate by default. The
        # same is also true for CUDA.
        config.check_home(conf, 'openmm')

        if openmm_type == 'reference':
            config.require_lib(conf, 'OpenMM_d')

            env.AppendUnique(CPPDEFINES = ['HAVE_OPENMM'])

        elif openmm_type == 'cuda':
            config.check_home(conf, 'cuda')

            config.require_lib(conf, 'OpenMM_d')
            config.require_lib(conf, 'OpenMMCuda_d')
            config.require_lib(conf, 'cudart')

            env.AppendUnique(CPPDEFINES = ['HAVE_OPENMM'])

    # Gromacs
    gromacs = env.get('gromacs', 0)
    if gromacs:
        config.check_home(conf, 'gromacs', '', '')

        config.require_lib(conf, 'md')
        config.require_lib(conf, 'gmx')

        config.require_header(conf, 'gromacs/txtdump.h')
        config.require_header(conf, 'gromacs/names.h')

        env.AppendUnique(CPPDEFINES = ['HAVE_GROMACS'])


def configure(conf):
    env = conf.env

    configure_deps(conf)

    # Libprotomol paths
    if home:
        env.Append(CPPPATH = [home + '/src'])
        env.Prepend(LIBPATH = [home])

    # Library name
    lib = 'protomol'
    if env.get('fah', 0): lib += '-fah'

    config.require_lib(conf, lib)
