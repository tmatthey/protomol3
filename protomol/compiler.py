def compiler_add_opts():
    opts.AddOptions(
        EnumOption('mode', 'Set build mode', 'debug',
                   allowed_values = ('debug', 'release')),
        BoolOption('optimize', 'Set to 1 to force optimizations', 1),
        BoolOption('gui', 'Set to 1 if using the GUI', 0),
        BoolOption('debug', 'Set to 1 to force debug options', 0),
        BoolOption('strict', 'Set to 0 to disable strict options', 1),
        BoolOption('depends', 'Set to 1 to output dependency files', 0),
        EnumOption('compiler', 'Select compiler', 'default',
                   allowed_values = ('default', 'intel', 'gnu', 'mingw', 'msvc',
                                     'linux-mingw', 'aix', 'posix', 'hp', 'sgi',
                                     'sun')))

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

    gui = int(env.get('gui', 0))    
    strict = int(env.get('strict', 1))
    depends = int(env.get('depends', 0))
    compiler = env.get('compiler', 0)

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
        env.Append(CCFLAGS = '-EHsc /Zp')


    # Debug flags
    if debug:
        if env['CC'] == 'cl':
            env.Append(CCFLAGS = '-W1 /MTd')
        elif env['CC'] == 'gcc':
            env.Append(CCFLAGS = '-g -Wall')
            if strict: env.Append(CCFLAGS = '-Werror')

        env.Append(CPPDEFINES = ['DEBUG'])

    else:
        if env['CC'] == 'gcc':
            env.Append(LINKFLAGS = '--strip-all')
        elif env['CC'] == 'cl':
            env.Append(CCFLAGS = '/MT')


    # Optimizations
    if optimize:
        if env['CC'] in ['icc', 'icpc']:
            env.Append(CCFLAGS = '-O -finline-functions')
        elif env['CC'] == 'gcc':
            env.Append(CCFLAGS =
                       '-O9 -ffast-math -finline-functions -funroll-loops')
            #env.Append(CCFLAGS = ['-msse2 -mfpmath=sse']);
        elif env['CC'] == 'cl':
            env.Append(CCFLAGS = '/Ox /GL')


    # GUI
    if gui:
        env.Append(CCFLAGS = '-DHAVE_GUI')

    # Dependency files
    if depends and env['CC'] == 'gcc':
        env.Append(CCFLAGS = '-MMD -MF ${TARGET}.d')


    # C mode
    if c99_mode:
        if env['CC'] == 'gcc':
            env.Append(CFLAGS = '-std=gnu99')
        elif env['CC'] == 'cl':
            env.Append(CFLAGS = '/TP') # C++ mode
