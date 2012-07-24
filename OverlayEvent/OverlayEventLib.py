#$Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['OverlayEvent'],
                 package = 'OverlayEvent')
    env.Tool('geometryLib')
    env.Tool('identsLib')
    env.Tool('addLibrary', library = env['clhepLibs'])
    env.Tool('addLibrary', library = env['gaudiLibs'])

def exists(env):
    return 1;
