# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['digiRootData'])
        if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'digiRootData') 
    env.Tool('commonRootDataLib')
    env.Tool('addLibrary', library = env['ldfLibs'])
    if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
        env.Tool('findPkgPath', package = 'enums') 

def exists(env):
    return 1;
