# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['digiRootData'])
    env.Tool('commonRootDataLib')
    env.Tool('addLibrary', library = env['ldfLibs'])
def exists(env):
    return 1;
