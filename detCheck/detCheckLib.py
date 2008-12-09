# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['detCheck'])
    env.Tool('addLibrary', library = env['xercesLibs'])
    env.Tool('detModelLib')
    env.Tool('identsLib')
def exists(env):
    return 1;
