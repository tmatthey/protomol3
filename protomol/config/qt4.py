import os
from SCons.Script import *

def configure(conf, libs = []):
    env = conf.env

    home = None
    if os.environ.has_key('QT4_HOME'):
        home = os.environ['QT4_HOME']
    elif env['PLATFORM'] == 'posix':
        home = '/usr'

    if home is None:
        raise Exception, 'Please set QT4_HOME'

    env.Append(LIBPATH = [home + '/lib'])

    if env['PLATFORM'] == 'win32':
        env.Append(CPPPATH = [home + '/include'])
    else:
        env.Append(CPPPATH = [home + '/include/qt4'])

    for lib in libs:
        if env['PLATFORM'] == 'win32':
            include = home + '/include/' + lib
        else:
            include = home + '/include/qt4/' + lib

        env.Append(CPPPATH = [include])

        if not conf.CheckCXXHeader(lib):
            raise Exception, 'QT4 library "' + lib + \
                '" header not found at ' + include

        if env['PLATFORM'] == 'win32': lib = lib + '4'

        if not conf.CheckLib(lib):
            raise Exception, 'QT4 library "' + lib + '" not found'
