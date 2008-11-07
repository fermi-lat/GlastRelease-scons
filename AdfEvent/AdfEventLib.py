# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['AdfEvent'])
    env.Tool('facilitiesLib')
    env.Tool('addLibrary', library = env['gaudiLibs'])
def exists(env):
    return 1;

