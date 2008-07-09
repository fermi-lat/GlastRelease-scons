# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['CalibData'])
    env.Tool('identsLib')

def exists(env):
    return 1;
