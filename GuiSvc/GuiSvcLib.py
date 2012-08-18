# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['GuiSvc'])
        if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'GuiSvc') 
    env.Tool('guiLib')
    env.Tool('addLibrary', library = env['gaudiLibs'])
def exists(env):
    return 1;
