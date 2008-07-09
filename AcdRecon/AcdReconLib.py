#$Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['AcdRecon'])
    env.Tool('EventLib')
    env.Tool('guiLib')
    env.Tool('identsLib')
    env.Tool('AcdUtilLib')

def exists(env):
    return 1;
