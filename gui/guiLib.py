# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['gui'])
    env.Tool('addLibrary', library = ['guisystem'])
    if env['PLATFORM'] != 'win32':
        env.Tool('addLibrary', library = env['x11Libs'])
    else:
        env.Tool('addLibrary', library = ['gdi32', 'winspool', 'comdlg32', 'shell32'])
    env.Tool('addLibrary', library = env['clhepLibs'])
def exists(env):
    return 1;
