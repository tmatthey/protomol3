import os
import tarfile
import re
import platform
import time
from SCons.Script import *


exclude_pats = [
    r'\.svn', r'\.sconsign.dblite', r'\.sconf_temp', r'.*~', r'.*\.o',
    r'.*\.obj'
    ]

distver = None
distrev = '-%(version)s-%(system)s-%(bits)s'

def add_vars(vars):
    vars.Add('dist_version', 'Set dist file version', None)
    vars.Add('dist_revision', 'Set dist file revision info', distrev)


def modify_targets(target, source, env):
    vars = {
        'machine' : platform.machine(),
        'bits' : platform.architecture()[0],
        'system' : platform.system(),
        'release' : platform.release(),
        'version' : env.get('dist_version', distver),
        'date' : time.strftime('%Y%m%d'),
        }

    rev = env.get('dist_revision', distrev)
    target = ((str(target[0]) + rev) % vars % env._dict) + '.tar.bz2'

    return [target, source]


def build_function(target, source, env):
    target = str(target[0])
    distname = os.path.splitext(os.path.splitext(target)[0])[0]

    tar = tarfile.open(target, mode = 'w:bz2')

    exclude_re = re.compile('^((' + ')|('.join(exclude_pats) + '))$')

    source = map(lambda x: str(x), source)

    def exclude(path):
        return exclude_re.match(os.path.split(path)[1]) != None

    def add_file(path):
        tar.add(path, distname + '/' + path, exclude = exclude)

    for src in source:
        if os.path.isdir(src):
            for root, dirs, names in os.walk(src):
                for dir in dirs:
                    if exclude(dir): dirs.remove(dir)

                for name in names:
                    add_file(os.path.join(root, name))
        else: add_file(src)

    return None


def configure(conf, default_distver, default_distrev = None,
              no_default_excludes = False, excludes = []):
    global exclude_pats, distver, distrev

    env = conf.env

    distver = default_distver
    if default_distrev: distrev = default_distrev

    if no_default_excludes: exclude_pats = []
    exclude_pats += excludes

    bld = Builder(action = build_function,
                  source_factory = SCons.Node.FS.Entry,
                  source_scanner = SCons.Defaults.DirScanner,
                  emitter = modify_targets)
    env.Append(BUILDERS = {'TarBZ2Dist' : bld})

    return True
