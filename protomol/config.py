# libFAH
if env['PLATFORM'] == 'win32':
    env.Append(LIBS = 'wsock32')
else:
    env.Append(LIBS = 'pthread')


# BOOST_HOME
if os.environ.has_key('BOOST_HOME'):
    env.Append(CPPPATH = [os.environ['BOOST_HOME']])


have_libfah=0
if os.environ.has_key('LIBFAH_HOME'):
    env.Append(CPPPATH = [os.environ['LIBFAH_HOME']])
    env.Append(LIBPATH = [os.environ['LIBFAH_HOME']])

if (conf.CheckLib('fah') and conf.CheckCXXHeader('fah/core/ChecksumDevice.h')):
    env.Append(CPPDEFINES = ['HAVE_LIBFAH'])
    have_libfah=1


# libbip2
have_bzip2=0
if os.environ.has_key('LIBBZ2_HOME'):
    env.Append(LIBPATH = [os.environ['LIBBZ2_HOME']])
    env.Append(CPPPATH = [os.environ['LIBBZ2_HOME'] + '/src'])

if not conf.CheckLib('bz2'):
    print 'libbzip2 not found.  Please set LIBBZ2_HOME'

else: have_bzip2 = 1


# boost::iostreams
if os.environ.has_key('BOOST_HOME'):
    env.Append(CPPPATH = [os.environ['BOOST_HOME']])

have_boost_iostreams=0
if (conf.CheckCXXHeader('boost/iostreams/stream.hpp')):
    have_boost_iostreams=1


# pthreads
if env['PLATFORM'] != 'win32' and have_libfah:
    if not conf.CheckLib('pthread'):
        print ('Need libpthreads for Folding@Home GUI server on non-win32 ' +
               'platforms')
        Exit(1)


if env['CC'] == 'gcc':
    if os.environ.has_key('ATLAS_HOME'):
        env.Append(LIBPATH = [os.environ['ATLAS_HOME']])
    else:
        env.Append(LIBPATH = ['/usr/lib/atlas'])


# LAPACK
if os.environ.has_key('LAPACK_HOME'):
    env.Append(CPPPATH = [os.environ['LAPACK_HOME']])
    env.Append(LIBPATH = [os.environ['LAPACK_HOME']])
  
have_lapack = 0
if conf.CheckLib('lapack'):
    env.Append(CPPDEFINES = ['HAVE_LAPACK'])
    have_lapack = 1


if env['CC'] == 'gcc' and have_lapack:
    # BLAS
    if os.environ.has_key('BLAS_HOME'):
        env.Append(LIBPATH = [os.environ['BLAS_HOME']])

    have_blas = 0  
    if conf.CheckLib('blas'):
        have_blas = 1


    # G2C
    if os.environ.has_key('G2C_HOME'):
        env.Append(LIBPATH = [os.environ['G2C_HOME']])

    have_g2c = 0  
    if conf.CheckLib('g2c'):
        have_g2c = 1
