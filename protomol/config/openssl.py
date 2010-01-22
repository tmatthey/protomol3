import os
from SCons.Script import *

def configure(conf):
    env = conf.env

    if os.environ.has_key('OPENSSL_HOME'):
        home = os.environ['OPENSSL_HOME']
        if os.path.exists(home + '/inc32') and os.path.exists(home + '/out32'):
            env.Append(CPPPATH = [home + '/inc32'])
            env.Append(LIBPATH = [home + '/out32'])
        else:
            env.Append(CPPPATH = [home + '/include'])
            env.Append(LIBPATH = [home + '/lib'])

    if (conf.CheckCHeader('openssl/ssl.h') and
        (conf.CheckLib('ssl') and conf.CheckLib('crypto')) or
        (conf.CheckLib('ssleay32') and conf.CheckLib('libeay32'))):

        if env['PLATFORM'] == 'posix':
            if not int(env.get('static', 0)): conf.CheckLib('dl')

        if env['PLATFORM'] == 'win32':
            for lib in ['wsock32', 'advapi32', 'gdi32', 'user32']:
                if not conf.CheckLib(lib):
                    raise Exception, 'openssl needs ' + lib

        env.Append(CPPDEFINES = ['HAVE_OPENSSL'])
        return True

    else: raise Exception, 'Need openssl'

