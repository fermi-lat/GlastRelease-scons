# $Header
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['AcdDigi'])
    env.Tool('AcdUtilLib')
    env.Tool('EventLib')
    env.Tool('GlastSvcLib')
def exists(env):
    return 1;
