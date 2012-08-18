# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['calibUtil'])
        if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'calibUtil') 
    env.Tool('facilitiesLib')
    env.Tool('xmlBaseLib')
    env.Tool('rdbModelLib')
    env.Tool('addLibrary', library = env['xercesLibs'])
def exists(env):
    return 1;
