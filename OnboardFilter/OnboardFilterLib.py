# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['OnboardFilter'])
 
    env.Tool('addLibrary', library = ['dl'])
    env.Tool('addLibrary', library = ['pthread'])
    env.Tool('addLibrary', library = env['obfLibs'])
    env.Tool('addLibrary', library = env['gaudiLibs'])
    env.Tool('OnboardFilterTdsLib')
    env.Tool('addLibrary', library = env['rootLibs'])
def exists(env):
    return 1;
