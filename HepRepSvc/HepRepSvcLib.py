# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['HepRepSvc'])
        if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'HepRepSvc') 
    env.Tool('CalUtilLib')
    env.Tool('CalibDataLib')
    env.Tool('GlastSvcLib')
    env.Tool('AcdUtilLib')
    env.Tool('addLibrary', library = env['clhepLibs'])
    if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
        env.Tool('findPkgPath', package = 'TkrUtil') 
        env.Tool('findPkgPath', package = 'FluxSvc') 
        env.Tool('findPkgPath', package = 'astro') 
        env.Tool('findPkgPath', package = 'facilities') 
        env.Tool('findPkgPath', package = 'ntupleWriterSvc') 
        env.Tool('findPkgPath', package = 'RootIo') 
        env.Tool('findPkgPath', package = 'rootUtil') 
        env.Tool('findPkgPath', package = 'OnboardFilterTds') 
def exists(env):
    return 1;
