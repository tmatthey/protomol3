import platform
import config

def MakePackages(env, target = None, source = None, progs = [], libs = [],
                 docs = [], screesavers = [], services = [], **kw):
    
    return []

    tool_deps = ['qt4','deb','packaging2']
    config.tools.load(env, tool_deps)
    targets = []

    name = kw['NAME'].lower()
    kw['NAME'] = name

    if platform.system() == 'Linux':
        source = []
        for x in progs:
            source.append(env.InstallAs('/usr/bin/' + str(x), str(x)))

        for x in docs:
            path = '/usr/share/docs/' + name + '/' + str(x)
            source.append(env.InstallAs(path, x))

        # Debian
        pkg = env.Package(PACKAGETYPE = 'deb', source = source, **kw)
        targets.append(pkg)

        # RPM
        #pkg = env.Package(PACKAGETYPE = 'rpm', RPMFLAGS = '-tb',
        #                  source = source, **kw)
        #targets.append(pkg)

    return targets


def configure(conf):
    conf.env.Append(BUILDERS = {'MakePackages' : MakePackages})

    return True
