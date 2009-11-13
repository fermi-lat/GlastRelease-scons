# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['TkrDigi'])
    env.Tool('facilitiesLib')
    env.Tool('EventLib')
    env.Tool('GlastSvcLib')
    env.Tool('TkrUtilLib')
    env.Tool('addLibrary', library = env['gaudiLibs'])
    env.Tool('addLibrary', library = env['clhepLibs'])
def exists(env):
    return 1;
