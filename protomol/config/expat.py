import os
from SCons.Script import *

def configure(conf):
    env = conf.env

    libname = 'expat'

    if os.environ.has_key('EXPAT_INCLUDE'):
        env.Append(CPPPATH = [os.environ['EXPAT_INCLUDE']])

    if os.environ.has_key('EXPAT_LIBPATH'):
        env.Append(LIBPATH = [os.environ['EXPAT_LIBPATH']])

    if env['PLATFORM'] == 'win32':
        libname = 'lib' + libname + 'MT'
        env.Append(CPPDEFINES = ['XML_STATIC'])

    if conf.CheckCHeader('expat.h') and conf.CheckLib(libname):
        env.Append(CPPDEFINES = ['HAVE_EXPAT'])

    else: raise Exception, 'Need expat'

    return True
