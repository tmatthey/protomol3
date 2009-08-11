execfile('configfuncs.py')

def add_fah_options():
    opts.AddOptions(
        BoolOption('fah', 'Set to 1 to build library for Folding@Home', 0)
    )

def fah_configure():
    # libfah
    libfah_home = check_envvar( 'LIBFAH_HOME', True )

    env.Append(CPPPATH = [libfah_home])
    env.Append(LIBPATH = [libfah_home])

    if check_library('fah', True) and check_header('fah/ID.h', True):
        env.Append(CPPDEFINES = ['HAVE_LIBFAH'])

    # boost
    boost_configure(
        conf, ['version', 'iostreams/stream'], ['iostreams', 'system', 'filesystem']
    )

    # boost::iostreams
    check_header( 'boost/iostreams/stream.hpp', True )

    # libbip2
    libbzip2_home = check_envvar( 'LIBBZIP2_HOME', True )

    env.Append(LIBPATH = [libbzip2_home])
    env.Append(CPPPATH = [libbzip2_home])

    check_library( 'bz2', True )

    # pthreads
    if env['PLATFORM'] != 'win32':
        check_library('pthread')

    env.Append(CPPDEFINES = ['BUILD_FOR_FAH', 'HAVE_LIBFAH'])
