# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['OnboardFilterTds'])
    env.Tool('EventLib')
    env.Tool('addLibrary', library=env['obfLibs'])
def exists(env):
    return 1;
