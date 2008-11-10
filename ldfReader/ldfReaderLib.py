# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = 'ldfReader')
    env.Tool('addLibrary', library = env['ldfLibs'])
    env.Tool('addLibrary', library = env['cfitsioLibs'])
    env.Tool('addLibrary', library = env['obfLibs'])
    env.Tool('facilitiesLib')
    env.Tool('astroLib')
    env.Tool('lsfDataLib')
def exists(env):
    return 1;
