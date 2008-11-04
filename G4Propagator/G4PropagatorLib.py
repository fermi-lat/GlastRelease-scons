# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['G4Propagator'])
    env.Tool('G4GeneratorLib')
def exists(env):
    return 1;
