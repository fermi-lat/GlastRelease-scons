#$Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['mcRootData'])
    env.Tool('commonRootDataLib')
def exists(env):
    return 1;
