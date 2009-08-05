# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['FluxSvc'])
    env.Tool('celestialSourcesLib')
    env.Tool('tipLib')
    env.Tool('fluxLib')
    env.Tool('TriggerLib')
    env.Tool('addLibrary', library = env['clhepLibs'])
def exists(env):
    return 1;
