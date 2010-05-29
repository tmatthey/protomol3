import os
import sys
import tarfile
import re
import platform
import time
from SCons.Script import *


exclude_pats = [
    r'\.svn', r'\.sconsign.dblite', r'\.sconf_temp', r'.*~', r'.*\.o',
    r'.*\.obj'
    ]

dist_ver = None
dist_build = '-%(version)s-%(system)s-%(bits)s'


def find_files(path, exclude = None):
    if not os.path.exists(path): return []
    if exclude:
        dir, filename = os.path.split(path)
        if exclude.match(filename): return []
    if not os.path.isdir(path): return [path]

    files = []
    for f in os.listdir(path):
        files += find_files(path + '/' + f, exclude)

    return files


def add_vars(vars):
    vars.Add('dist_version', 'Set dist file version', None)
    vars.Add('dist_build', 'Set dist file build info', dist_build)


def modify_targets(target, source, env):
    vars = {
        'machine' : platform.machine(),
        'bits' : platform.architecture()[0],
        'system' : platform.system(),
        'release' : platform.release(),
        'version' : env.get('dist_version', dist_ver),
        'date' : time.strftime('%Y%m%d'),
        }

    build = env.get('dist_build', dist_build)
    if len(build) and build[0] != '-': build = '-' + build

    if len(dist_ver): build = '-' + dist_ver + build

    target = ((str(target[0]) + build) % vars % env._dict) + '.tar.bz2'

    # Write 'dist.txt'
    f = None
    try:
        f = open('dist.txt', 'w')
        f.write(target)
    finally:
        if f is not None: f.close()

    return [target, source]


def build_function(target, source, env):
    target = str(target[0])
    dist_name = os.path.splitext(os.path.splitext(target)[0])[0]

    tar = tarfile.open(target, mode = 'w:bz2')

    exclude_re = re.compile('^((' + ')|('.join(exclude_pats) + '))$')

    source = map(lambda x: str(x), source)

    for src in source:
        for file in find_files(src, exclude_re):
            tar.add(file, dist_name + '/' + file)

    return None


def configure(conf, version = '', no_default_excludes = False, excludes = []):
    global exclude_pats, dist_ver, dist_build

    env = conf.env

    dist_ver = version

    if no_default_excludes: exclude_pats = []
    exclude_pats += excludes

    bld = Builder(action = build_function,
                  source_factory = SCons.Node.FS.Entry,
                  source_scanner = SCons.Defaults.DirScanner,
                  emitter = modify_targets)
    env.Append(BUILDERS = {'TarBZ2Dist' : bld})

    return True
