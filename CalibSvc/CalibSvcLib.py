# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['CalibSvc'])
    env.Tool('CalibDataLib')
    env.Tool('xmlBaseLib')
    env.Tool('calibUtilLib')
    env.Tool('EventLib')
    env.Tool('commonRootDataLib')
    env.Tool('calibRootDataLib')
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['rootGuiLibs'])
    env.Tool('LdfEventLib')
def exists(env):
    return  1;
