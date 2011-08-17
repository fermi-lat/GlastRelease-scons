# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['G4Generator'])
    env.Tool('fluxLib')
    env.Tool('EventLib')
    env.Tool('addLibrary', library = env['geant4Libs'])
    env.Tool('GlastSvcLib')
    env.Tool('geometryLib')
    env.Tool('addLibrary', library = env['clhepLibs'])
def exists(env):
    return 1;
