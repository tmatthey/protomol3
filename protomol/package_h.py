REV="unknown"
REPO="unknown"

def c_str_escape(str):
   return str.replace('\\', '\\\\').replace('"', '\\"').replace('\n', '\\n')

if os.path.exists( '.svn/entries' ):
  f = open('.svn/entries')
  try:
      i = 0
      for line in f:
          l = c_str_escape(line.strip())
          if l.startswith("http"):
              REPO = line.strip()
              break
          else:
              REV = line.strip()
  finally:
      f.close()

tags = [
    ['VERSION', '3.3'],
    ['BUGREPORT', '\nprotomol@cse.nd.edu, joseph@cauldrondevelopment.com'],
    ['HOMEPAGE', 'http://protomol.sourceforge.net/'],
    ['CITE',
     'T. Matthey, T. Cickovski, S. S. Hampton, A. Ko, Q. Ma,\n' +
     'M. Nyerges, T. Raeder, T. Slabach, and J. A. Izaguirre.\nProtoMol: ' +
     'An object-oriented framework for prototyping novel algorithms\nfor ' +
     'molecular dynamics.\nACM Trans. Math. Softw., 30(3):237 265, 2004.'],
]

computed_tags = [
    ['REVISION', REV],
    ['SOURCE_REPO', '\n' + REPO],
    ['COMPILER', env.subst('$CXX')],
    ['COMPILER_VERSION', env.subst('$CXXVERSION')],
    ['COMPILER_FLAGS', env.subst('$CXXFLAGS $_CPPDEFFLAGS')],
    ['COMPILER_LIBS', env.subst('$LIBS')],
    ['BUILT_BY', os.environ.get('USER', '')],
    ['PLATFORM', env.subst('$PLATFORM')],
]

f = open('src/protomol/package.h', 'w')
try:
    f.write('#ifndef PROTOMOL_PACKAGE_H\n')
    f.write('#define PROTOMOL_PACKAGE_H\n\n')
    f.write('/************************************************************\n')
    f.write(' * NOTE: This is an autogenerated file.  Please do not edit *\n')
    f.write(' *       directly or check in to revision control.          *\n')
    f.write(' ************************************************************/\n')
    f.write('\n')

    for tag in tags + computed_tags:
        f.write('#define PACKAGE_' + tag[0] + ' ')
        if tag[1].find('\n') != -1:
            i = 0
            for l in tag[1].split('\n'):
                if l != '':
                    i = i + 1
                    if i > 1: f.write('\\n"')
                    f.write('\\\n  "' + l)
            f.write('"\n')
        else:
            f.write('"' + tag[1] + '"\n')

    f.write('\n#endif // PROTOMOL_PACKAGE_H\n')
    env.Append(CPPDEFINES = ['HAVE_PACKAGE_H'])

finally:
    f.close()
