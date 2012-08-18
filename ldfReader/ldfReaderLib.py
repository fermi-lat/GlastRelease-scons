# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['ldfReader'])
        if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'ldfReader') 
    env.Tool('addLibrary', library = env['ldfLibs'])
    env.Tool('facilitiesLib')
    env.Tool('eventFileLib')
    env.Tool('lsfDataLib')
    env.Tool('astroLib')
    if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
        env.Tool('findPkgPath', package = 'enums') 

def exists(env):
    return 1;
