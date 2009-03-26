def compiler_add_opts():
    opts.AddOptions(
        EnumOption('mode', 'Set build mode', 'debug',
                   allowed_values = ('debug', 'release')),
        BoolOption('optimize', 'Set to 1 to force optimizations', 0),
        BoolOption('debug', 'Set to 1 to force debug options', 0),
        BoolOption('strict', 'Set to 0 to disable strict options', 1),
        BoolOption('threaded', 'Set to 1 to enable thread support', 1),
        BoolOption('profile', 'Set to 1 to enable profiler', 0),
        BoolOption('depends', 'Set to 1 to output dependency files', 0),
        BoolOption('distcc', 'Set to 1 to enable distributed builds', 0),
        BoolOption('ccache', 'Set to 1 to enable cached builds', 0),
        EnumOption('platform', 'Override default platform', '',
                   allowed_values = ('', 'win32', 'posix')),
        EnumOption('cxxstd', 'Set C++ language standard', 'c++98',
                   allowed_values = ('c++98', 'c++0x')),
        EnumOption('compiler', 'Select compiler', 'default',
                   allowed_values = ('default', 'gnu', 'intel', 'mingw', 'msvc',
                                     'linux-mingw', 'aix', 'posix', 'hp', 'sgi',
                                     'sun')))

    opts.Add('static', 'Link to static libraries', 0)
    opts.Add('num_jobs', 'Set the concurrency level.', -1)


def compiler_configure(c99_mode = 1):
    global depends, compiler, strict, optimize, debug, mode

    if env.GetOption('clean'):
        return

    # Get options
    mode = env.get('mode', 'debug')

    if env.has_key('optimize'): optimize = int(env['optimize'])
    else: optimize = mode == 'release'
    
    if env.has_key('debug'): debug = int(env['debug'])
    else: debug = mode == 'debug'
    
    strict = int(env.get('strict', 1))
    threaded = int(env.get('threaded', 1))
    profile = int(env.get('profile', 0))
    depends = int(env.get('depends', 0))
    compiler = env.get('compiler', 0)
    distcc = env.get('distcc', 0)
    ccache = env.get('ccache', 0)
    cxxstd = env.get('cxxstd', 'c++0x')
    platform = env.get('platform', '')
    static = int(env.get('static', 0))
    num_jobs = env.get('num_jobs', -1)

    if platform != '':
        env.Replace(PLATFORM = platform)

    # Select compiler
    if compiler:
        if compiler == 'gnu':
            Tool('gcc')(env)
            Tool('g++')(env)

        elif compiler == 'intel':
            Tool('intelc')(env)
            env['ENV']['INTEL_LICENSE_FILE'] = (
                os.environ.get("INTEL_LICENSE_FILE", ''))

        elif compiler == 'linux-mingw':
            env.Replace(CC = 'i586-mingw32msvc-gcc')
            env.Replace(CXX = 'i586-mingw32msvc-g++')
            env.Replace(RANLIB = 'i586-mingw32msvc-ranlib')
            env.Replace(PROGSUFFIX = '.exe')

        elif compiler == 'posix':
            Tool('cc')(env)
            Tool('cxx')(env)
            Tool('link')(env)
            Tool('ar')(env)
            Tool('as')(env)
            
        elif compiler in Split('hp sgi sun aix'):
            Tool(compiler + 'cc')(env)
            Tool(compiler + 'c++')(env)
            Tool(compiler + 'link')(env)
            
            if compiler in Split('sgi sun'):
                Tool(compiler + 'ar')(env)

        elif compiler != 'default':
            Tool(compiler)(env)


    print "Compiler: " + env['CC']
    print "Platform: " + env['PLATFORM']


    # Options
    if env['CC'] == 'cl':
        env.Append(CCFLAGS = ['/EHsc', '/Zp'])


    # Profiler flags
    if profile:
        if env['CC'] == 'gcc':
            env.Append(CCFLAGS = ['-pg'])
            env.Append(LINKFLAGS = ['-pg'])


    # Debug flags
    if debug:
        if env['CC'] == 'cl':
            env.Append(CCFLAGS = ['/W1'])
            env.Append(LINKFLAGS = ['/DEBUG', '/MAP:${TARGET}.map'])
            env['PDB'] = '${TARGET}.pdb'

        elif env['CC'] == 'gcc':
            env.Append(CCFLAGS = ['-ggdb', '-Wall'])
            if strict: env.Append(CCFLAGS = ['-Werror'])
            env.Append(LINKFLAGS = ['-rdynamic']) # for backtrace

        env.Append(CPPDEFINES = ['DEBUG'])

    else:
        if env['CC'] == 'gcc':
            env.Append(LINKFLAGS = ['--strip-all'])


    # Optimizations
    if optimize:
        if env['CC'] in ['icc', 'icpc']:
            env.Append(CCFLAGS = ['-O', '-finline-functions', '-funroll-loops'])
        elif env['CC'] == 'gcc':
            env.Append(CCFLAGS =
                       ['-O9', '-ffast-math', '-funroll-loops'])
            #env.Append(CCFLAGS = ['-msse2 -mfpmath=sse']);
        elif env['CC'] == 'cl':
            env.Append(CCFLAGS = ['/Ox', '/GL'])
            env.Append(LINKFLAGS = ['/LTCG'])
            env.Append(ARFLAGS = ['/LTCG'])


    # Dependency files
    if depends and env['CC'] == 'gcc':
        env.Append(CCFLAGS = ['-MMD -MF ${TARGET}.d'])


    # C mode
    if c99_mode:
        if env['CC'] == 'gcc':
            env.Append(CFLAGS = ['-std=c99'])
            env.Append(CXXFLAGS = ['-std=' + cxxstd])
        elif env['CC'] == 'cl':
            env.Append(CFLAGS = ['/TP']) # C++ mode


    # Threads
    if threaded:
        if env['CC'] == 'gcc':
            if not conf.CheckLib('pthread'):
                print 'Need pthreads'
                Exit(1)

            env.Append(LINKFLAGS = ['-pthread'])
            env.Append(CPPDEFINES = ['_REENTRANT'])

        elif env['CC'] == 'cl':
            if debug:
                env.Append(CCFLAGS = ['/MTd'])
            else:
                env.Append(CCFLAGS = ['/MT'])


    # static
    if static:
        if env['CC'] == 'gcc':
            env.Append(LINKFLAGS = ['-static'])


    # Num jobs default
    default_num_jobs = 1

    # distcc
    if distcc:
        default_num_jobs = 2
        env.Replace(CC = 'distcc ' + env['CC'])
        env.Replace(CXX = 'distcc ' + env['CXX'])

    # cccache
    if ccache:
        env.Replace(CC = 'ccache ' + env['CC'])
        env.Replace(CXX = 'ccache ' + env['CXX'])


    # Num jobs
    if num_jobs == -1:
        if os.environ.has_key('SCONS_JOBS'):
            num_jobs = int(os.environ.get('SCONS_JOBS', num_jobs))
        else:
            num_jobs = default_num_jobs

    SetOption('num_jobs', num_jobs)
    print "running with -j", GetOption('num_jobs')

