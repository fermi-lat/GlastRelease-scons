#$Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['CalUtil'])
        if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'CalUtil') 
    env.Tool('identsLib')
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['gaudiLibs'])
    env.Tool('xmlBaseLib')    
    env.Tool('EventLib')
    if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
        env.Tool('findPkgPath', package = 'GlastSvc') 
def exists(env):
    return 1;
