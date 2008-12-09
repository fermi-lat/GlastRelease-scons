# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['LdfConverter'])
    env.Tool('ldfReaderLib')	
    env.Tool('GlastSvcLib')
    env.Tool('EventLib')
def exists(env):
    return 1;
