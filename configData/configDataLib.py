# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['configData'])
    env.Tool('mootCoreLib')
    env.Tool('commonRootDataLib')
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['pythonLibs'])
    env.Tool('lsfDataLib')
def exists(env):
    return 1;
