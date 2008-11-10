# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['Interleave'])
    env.Tool('addLibrary', library = env['clhepLibs'])
    env.Tool('FluxSvcLib')
    env.Tool('CalUtilLib')
    env.Tool('CalibDataLib')
    env.Tool('GlastSvcLib')	
    env.Tool('EventLib')
def exists(env):
    return 1;
