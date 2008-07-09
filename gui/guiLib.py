# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['gui', 'guisystem'])
        env.Tool('addLibrary', library = env['x11Libs'])

def exists(env):
    return 1;
