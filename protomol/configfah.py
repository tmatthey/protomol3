# BOOST_HOME
if os.environ.has_key('BOOST_HOME'):
    env.Append(CPPPATH = [os.environ['BOOST_HOME']])

boost_configure(conf, ['version', 'iostreams/stream'], ['iostreams'])

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
    Exit(1)

else: have_bzip2 = 1


# boost::iostreams
have_boost_iostreams=0
if (conf.CheckCXXHeader('boost/iostreams/stream.hpp')):
    have_boost_iostreams=1


# pthreads
if env['PLATFORM'] != 'win32' and have_libfah:
    if not conf.CheckLib('pthread'):
        print ('Need libpthreads for Folding@Home GUI server on non-win32 ' +
               'platforms')
        Exit(1)
