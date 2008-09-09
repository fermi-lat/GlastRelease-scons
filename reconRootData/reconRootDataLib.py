# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['reconRootData'])
    env.Tool('commonRootDataLib')
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('facilitiesLib')
def exists(env):
    return 1;
