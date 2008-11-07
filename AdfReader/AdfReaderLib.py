#$Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['AdfReader'])
    env.Tool('facilitiesLib')
    env.Tool('AncillaryDataEventLib')
    env.Tool('AncillaryDataUtilLib')
    env.Tool('AdfEventLib')
def exists(env):
    return 1;
