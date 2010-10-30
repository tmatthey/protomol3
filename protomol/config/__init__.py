import sys
import os
import imp
from SCons.Script import *

modules = ['compiler', 'boost', 'python', 'xml', 'pthreads', 'expat', 'openssl',
           'qt4', 'dist', 'resources', 'build_info', 'packager', 'tools',
           'valgrind']

# Child modules
for name in modules:
    try:
        name = __name__ + '.' + name
        mod = __import__(name)
        sys.modules[__name__].__dict__[name] = mod
    except ImportError: pass


def require_lib(conf, lib):
    if not conf.CheckLib(lib):
        print 'Need library ' + lib
        Exit(1)


def load_conf_module(name, home = None):
    modname = 'config.' + name;

    if not modname in sys.modules:
        env_name = name.upper().replace('-', '_')
        env_home = env_name + '_HOME'

        if os.environ.has_key(env_home):
            home = os.environ[env_home]

        if home and not os.path.isdir(home):
            raise Exception, '$%s=%s is not a directory' % (env_home, home)

        filename = name + '-config.py'
        if home: filename = home + '/' + filename

        if os.path.exists(filename):
            try:
                imp.load_source(modname, filename)

                if home:
                    mod = sys.modules[modname]
                    mod.__dict__['home'] = home

            except Exception, e:
                raise Exception, "Error loading configuration file '%s' " \
                    "for '%s':\n%s" % (filename, name, e)
        else:
            raise Exception, "Configuration file '%s' not found " \
                "for '%s' not found please set %s" % (filename, name, env_home)

    return sys.modules[modname]



def call_single(mod, func, required = False, args = [], kwargs = {}):
    if type(mod) == str:
        try:
            mod = load_conf_module(mod)
        except:
            if required: raise
            return False

    try:
        return mod.__dict__[func](*args, **kwargs)

    except KeyError:
        if required:
            print "function '%s' not found in '%s'" % (func, mod.__name__)
            Exit(1)

        return False

    except Exception, e:
        if required:
            print e
            Exit(1)

        return False


def call_conf_module(mods, func, required = False, args = [], kwargs = {}):
    ret = True

    if not isinstance(mods, list): mods = [mods]
    for mod in mods:
        if not call_single(mod, func, required, args, kwargs):
            ret = False

    return ret


def get_deps(mods):
    if not isinstance(mods, list): mods = [mods]

    deps = set()
    for mod in mods:
        if not mod in deps:
            deps.add(mod)
            
            try:
                mod = load_conf_module(mod)
                mods += mod.__dict__['deps']

            except KeyError:
                pass

            except Exception:
                pass

    return list(deps)


def add_vars(mods, vars, required = False):
    if not isinstance(mods, list): mods = [mods]
    mods = get_deps(mods)

    return call_conf_module(mods, 'add_vars', required, [vars])


def configure(mods, conf, required = True, **kwargs):
    return call_conf_module(mods, 'configure', required, [conf], kwargs)


def make_env(deps, extra_vars = None, tool_deps=[]):
    configs = []
    
    if os.environ.has_key('SCONS_OPTIONS'):
        options = os.environ['SCONS_OPTIONS']
        if not os.path.exists(options):
            print 'options file "%s" set in SCONS_OPTIONS does not exist' % \
                options
            Exit(1)

        configs.append(options)

    if os.path.exists('default-options.py'):
        configs.append('default-options.py')
    if os.path.exists('options.py'): configs.append('options.py')

    vars = Variables(configs)
    add_vars(deps, vars)
    if extra_vars:
        for var in extra_vars: vars.AddVariables(var)
    env = Environment(variables = vars, ENV = os.environ)
    Help(vars.GenerateHelpText(env))

    try:
        tools.load(env, tool_deps)
    except NameError: pass
        
    Export('env')

    return env
