# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['gui'])
    env.Tool('addLibrary', library = ['guisystem'])
    env.Tool('addLibrary', library = env['x11Libs'])
    env.Tool('addLibrary', library = env['clhepLibs'])
def exists(env):
    return 1;
