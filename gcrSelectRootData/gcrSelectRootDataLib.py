# $Header$ 

def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['gcrSelectRootData'])
        if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'gcrSelectRootData') 
    env.Tool('commonRootDataLib')
    env.Tool('addLibrary', library = env['rootLibs'])
def exists(env):
    return  1;
