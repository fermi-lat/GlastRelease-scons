# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['AcdUtil', 'AcdUtilCommon'])
    env.Tool('CalibDataLib')
    env.Tool('geometryLib')
    env.Tool('xmlBaseLib')

def exists(env):
    return 1;
