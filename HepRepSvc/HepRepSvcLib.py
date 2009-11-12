# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['HepRepSvc'])
    env.Tool('CalUtilLib')
    env.Tool('CalibDataLib')
    env.Tool('GlastSvcLib')
    env.Tool('AcdUtilLib')
    env.Tool('addLibrary', library = env['clhepLibs'])
def exists(env):
    return 1;
