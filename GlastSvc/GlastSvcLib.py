#$Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['GlastSvc'])
    env.Tool('EventLib')
    env.Tool('xmlUtilLib')
    env.Tool('detModelLib')
    env.Tool('identsLib')
    env.Tool('facilitiesLib')

def exists(env):
    return 1;
