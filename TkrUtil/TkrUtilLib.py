# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['TkrUtil'])
        if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'TkrUtil') 
    env.Tool('CalibDataLib')
    env.Tool('xmlBaseLib')
    env.Tool('geometryLib')
    env.Tool('EventLib')
    if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
        env.Tool('findPkgPath', package='GlastSvc')
        env.Tool('findPkgPath', package='CalibSvc')
        env.Tool('findPkgPath', package='LdfEvent')

def exists(env):
    return 1;
