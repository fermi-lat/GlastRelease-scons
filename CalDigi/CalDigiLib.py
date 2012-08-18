#$Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['CalDigi'])
    env.Tool('configDataLib')
    env.Tool('EventLib')
    env.Tool('CalUtilLib')
    env.Tool('xmlBaseLib')
    env.Tool('addLibrary', library = env['gaudiLibs'])
    if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
        env.Tool('findPkgPath', package = 'CalXtalResponse') 
        env.Tool('findPkgPath', package = 'LdfEvent') 
        env.Tool('findPkgPath', package = 'Event') 
        env.Tool('findPkgPath', package = 'enums') 
	env.Tool('findPkgPath', package = 'ConfigSvc')
def exists(env):
    return 1;
