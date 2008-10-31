# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['FluxSvc'])
    env.Tool('celestialSourcesLib')
def exists(env):
    return 1;
