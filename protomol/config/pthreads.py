import os
from SCons.Script import *

def configure(conf):
    env = conf.env

    if os.environ.has_key('PTHREADS_HOME'):
        home = os.environ['PTHREADS_HOME']
        env.Append(CPPPATH = [home])
        env.Append(LIBPATH = [home])

    if conf.CheckCHeader('pthread.h') and conf.CheckLib('pthread'):
        env.Append(CPPDEFINES = ['HAVE_PTHREADS'])
    else: raise Exception, 'Need pthreads'
