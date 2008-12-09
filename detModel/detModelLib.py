# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['detModel'])
    env.Tool('addLibrary', library = env['clhepLibs'])
    env.Tool('xmlUtilLib')
def exists(env):
    return 1;
