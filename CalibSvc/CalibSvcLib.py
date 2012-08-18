# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['CalibSvc'])
        if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'CalibSvc') 
    env.Tool('CalibDataLib')
    env.Tool('xmlBaseLib')
    env.Tool('calibUtilLib')
    env.Tool('EventLib')
    env.Tool('commonRootDataLib')
    env.Tool('calibRootDataLib')
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['rootGuiLibs'])
    env.Tool('LdfEventLib')
    if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
        env.Tool('findPkgPath', package = 'MootSvc') 
def exists(env):
    return  1;
