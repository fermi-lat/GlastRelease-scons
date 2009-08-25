#$Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['CalDigi'])
    env.Tool('configDataLib')
    env.Tool('EventLib')
    env.Tool('CalUtilLib')
    env.Tool('xmlBaseLib')
    env.Tool('addLibrary', library = env['gaudiLibs'])
def exists(env):
    return 1;
