import os
import tarfile
import re
import platform
import time
import textwrap
from SCons.Script import *


exclude_pats = [r'\.svn', r'.*~']
exclude = None
ns = None
next_id = 0
col = 0


def write_string(output, s, newline = 0):
    global col

    l = len(s)

    if newline or col + l > 80:
        output.write('\n')
        col = 0

    i = s.rfind('\n')
    if i != -1: col = l - (i + 1)
    else: col += l

    output.write(s)


def write_resource(output, path, children = None, exclude = None):
    global next_id

    name = os.path.basename(path)
    if exclude != None and exclude.search(path) != None:
        return

    is_dir = os.path.isdir(path)
    id = next_id
    next_id += 1
    length = 0

    if is_dir:
        typeStr = 'Directory'
        child_resources = []
        
        for filename in os.listdir(path):
            write_resource(output, os.path.join(path, filename),
                           child_resources, exclude)

        length = len(child_resources)

        if length:
            write_string(output, 'const Resource *children%d[] = {' % id)
            for res in child_resources:
                write_string(output, '&resource%d,' % res)

            write_string(output, '0};\n')

    else:
        typeStr = 'File'
        f = open(path, 'r')

        write_string(output, 'const unsigned char data%d[] = {' % id)

        while True:
            buf = f.read(1024)

            length += len(buf)
            if len(buf) == 0: break

            for c in buf: write_string(output, '%d,' % ord(c))

            if len(buf) < 1024: break

        write_string(output, '0};\n')

    if children != None: children.append(id)

    output.write('extern const %sResource resource%d("%s", ' %
                 (typeStr, id, name))

    if is_dir: output.write('children%d' % id)
    else: output.write('(const char *)data%d, %d' % (id, length))

    output.write(');\n')


def build_function(target, source, env):
    global exclude, next_id, col, ns
    next_id = col = 0

    target = str(target[0])
    f = open(target, 'w')

    note = ('WARNING: This file was auto generated.  Please do NOT '
            'edit directly or check in to source control.')

    f.write(
        '/' + ('*' * 75) + '\\\n   ' +
        '\n   '.join(textwrap.wrap(note)) + '\n' +
        '\\' + ('*' * 75) + '/\n'
        '\n'
        '#include <fah/util/Resource.h>\n\n'
        )

    if ns:
        for namespace in ns.split('::'):
            f.write('namespace %s {\n' % namespace)

    for src in source:
        write_resource(f, str(src), exclude = exclude)
    
    if ns:
        for namespace in ns.split('::'):
            f.write('} // namespace %s\n' % namespace)

    f.close()

    return None


def configure(conf, namespace, no_default_excludes = False, excludes = []):
    global exclude, exclude_pats, ns

    env = conf.env

    ns = namespace

    if not no_default_excludes: excludes = list(excludes) + exclude_pats

    pattern = None
    for ex in excludes:
        if pattern == None: pattern = ''
        else: pattern += '|'
        pattern += '(%s)' % ex

    exclude = re.compile(pattern)

    bld = Builder(action = build_function,
                  source_factory = SCons.Node.FS.Entry,
                  source_scanner = SCons.Defaults.DirScanner)
    env.Append(BUILDERS = {'Resources' : bld})

    return True
