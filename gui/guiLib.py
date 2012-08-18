# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['gui'])
        env.Tool('addLibrary', library = ['guisystem'])
        if env['PLATFORM'] == "win32" and env.get('CONTAINERNAME','') == 'GlastRelease':
            env.Tool('findPkgPath', package = 'gui') 
    if env['PLATFORM'] != 'win32':
        env.Tool('addLibrary', library = env['x11Libs'])
    else:
        env.Tool('addLibrary', library = ['gdi32', 'winspool', 'comdlg32', 'shell32'])
    env.Tool('addLibrary', library = env['clhepLibs'])
    # For the time being don't need incsOnly branch since gui source
    # doesn't reference other packages
    #if kw.get('incsOnly', 0) == 1: 
    #    env.Tool('findPkgPath', package = ...

def exists(env):
    return 1;
