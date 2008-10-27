# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['digiRootData'])
    env.Tool('commonRootDataLib')
    env.Tool('facilitiesLib')
    env.Tool('addLibrary', library = env['ldfLibs'])
    env.Tool('lsfDataLib')
    env.Tool('addLibrary', library = env['obfLibs'])
    env.Tool('addLibrary', library = env['rootLibs'])
def exists(env):
    return 1;
