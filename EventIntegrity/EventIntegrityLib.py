# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['EventIntegrity'])
    env.Tool('GlastSvcLib')
    env.Tool('LdfEventLib')
def exists(env):
    return 1;
