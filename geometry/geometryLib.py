# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['geometry'])
        if env['PLATFORM'] == "win32" and env.get('CONTAINERNAME','') == 'GlastRelease':
            env.Tool('findPkgPath', package = 'geometry') 
    env.Tool('addLibrary', library = env['clhepLibs'])
    # For the time being don't need incsOnly branch since geometry source
    # doesn't reference other packages
    #if kw.get('incsOnly', 0) == 1: 
    #    env.Tool('findPkgPath', package = ...

def exists(env):
    return 1;
