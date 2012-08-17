#$Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['ConfigSvc'])
        if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'ConfigSvc') 

    env.Tool('configDataLib')
    env.Tool('MootSvcLib')
    env.Tool('EventLib')
    env.Tool('LdfEventLib')
    env.Tool('addLibrary', library = env['gaudiLibs'])
    if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
        env.Tool('findPkgPath', package = 'enums') 

def exists(env):
    return 1;
