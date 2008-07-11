# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['TkrUtil'])
    env.Tool('CalibDataLib')
    env.Tool('xmlBaseLib')
    env.Tool('geometryLib')
    env.Tool('EventLib')

def exists(env):
    return 1;
