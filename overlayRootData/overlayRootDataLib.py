# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['overlayRootData'])
    env.Tool('commonRootDataLib')
    env.Tool('facilitiesLib')
def exists(env):
    return 1;
