#$Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['DetDisplay'])
    env.Tool('TkrUtilLib')
    env.Tool('GlastSvcLib')
    env.Tool('geomrepLib')
def exists(env):
    return 1;
