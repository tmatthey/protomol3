import textwrap
from platform import *
from SCons.Script import *

ns = ''

def build_function(target, source, env):
    entries = '.svn/entries'
    if os.path.exists(entries):
        f = open(entries)
        for i in range(4):
            line = f.readline()

        revision = line.strip()
    else: revision = 'Unknown'


    target = str(target[0])
    f = open(target, 'w')

    note = ('WARNING: This file was auto generated.  Please do NOT '
            'edit directly or check in to source control.')

    f.write(
        '/' + ('*' * 75) + '\\\n   ' +
        '\n   '.join(textwrap.wrap(note)) + '\n' +
        '\\' + ('*' * 75) + '/\n'
        '\n'
        '#include <fah/Info.h>\n'
        '#include <fah/String.h>\n'
        '#include <fah/util/CompilerInfo.h>\n\n'
        )

    if ns:
        for namespace in ns.split('::'):
            f.write('namespace %s {\n' % namespace)

    if env.get('debug', False): mode = 'Debug'
    else: mode = 'Release'

    f.write(
        '  void addBuildInfo(const char *category) {\n'
        '    Info &info = Info::instance();\n'
        '\n'
        '    info.add(category, "Date", __DATE__);\n'
        '    info.add(category, "Time", __TIME__);\n'
        '    info.add(category, "Revision", "%s");\n'
        '    info.add(category, "Compiler", COMPILER);\n'
        '    info.add(category, "Options", "%s");\n'
        '    info.add(category, "Platform", "%s");\n'
        '    info.add(category, "Bits", String(COMPILER_BITS));\n'
        '    info.add(category, "Mode", "%s");\n'
        '  }\n' % (
            revision,
            ' '.join(env['CXXFLAGS']) + ' ' + ' '.join(env['CCFLAGS']),
            system() + ' ' + release(), mode,
            ))

    if ns:
        for namespace in ns.split('::'):
            f.write('} // namespace %s\n' % namespace)

    f.close()

    return None


def configure(conf, namespace):
    global ns

    ns = namespace

    env = conf.env

    bld = Builder(action = build_function,
                  source_factory = SCons.Node.FS.Entry,
                  source_scanner = SCons.Defaults.DirScanner)
    env.Append(BUILDERS = {'BuildInfo' : bld})

    return True
