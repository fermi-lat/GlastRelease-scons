# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['MootSvc'])
        if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'MootSvc') 

    env.Tool('addLibrary', library = env['gaudiLibs'])
    env.Tool('facilitiesLib')
    env.Tool('xmlBaseLib')
    env.Tool('EventLib')
    env.Tool('LdfEventLib')
    env.Tool('mootCoreLib')
    env.Tool('CalibDataLib')
def exists(env):
    return 1;
