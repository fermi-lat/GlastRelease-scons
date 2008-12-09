# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['Gleam'])
    env.Tool('OnboardFilterLib')
def exists(env):
    return 1;
