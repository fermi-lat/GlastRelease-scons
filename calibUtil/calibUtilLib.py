# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['calibUtil'])
    env.Tool('facilitiesLib')
    env.Tool('xmlBaseLib')
    env.Tool('rdbModelLib')
    env.Tool('addLibrary', library = env['xercesLibs'])
def exists(env):
    return 1;
