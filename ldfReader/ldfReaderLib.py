# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['ldfReader'])
    env.Tool('addLibrary', library = env['ldfLibs'])
    env.Tool('facilitiesLib')
    env.Tool('eventFileLib')
    env.Tool('lsfDataLib')
    env.Tool('astroLib')
def exists(env):
    return 1;
