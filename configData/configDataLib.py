# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['configData'])
    env.Tool('mootCoreLib')
    env.Tool('commonRootDataLib')
def exists(env):
    return 1;
