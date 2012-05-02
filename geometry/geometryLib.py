# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['geometry'])
    env.Tool('addLibrary', library = env['clhepLibs'])
def exists(env):
    return 1;
