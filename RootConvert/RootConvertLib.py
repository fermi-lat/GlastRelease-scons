# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['RootConvert'])
        if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'RootConvert') 
    env.Tool('EventLib')
    env.Tool('LdfEventLib')
    env.Tool('OnboardFilterTdsLib')
    env.Tool('addLibrary', library = env['obfLibs'])
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['gaudiLibs'])
    env.Tool('reconRootDataLib')
    env.Tool('digiRootDataLib')
    env.Tool('mcRootDataLib')
    env.Tool('commonRootDataLib')
    env.Tool('gcrSelectRootDataLib')
    env.Tool('AncillaryDataEventLib')
def exists(env):
    return 1;
