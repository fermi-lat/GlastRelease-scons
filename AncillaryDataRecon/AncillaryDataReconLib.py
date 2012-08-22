# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['AncillaryDataRecon'])
        if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'AncillaryDataRecon') 
    env.Tool('facilitiesLib')
    env.Tool('CalibDataLib')
    env.Tool('AncillaryDataEventLib')
    env.Tool('AncillaryDataUtilLib')
    env.Tool('AdfReaderLib')
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['rootGuiLibs'])
    env.Tool('addLibrary', library = env['clhepLibs'])
    if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
        env.Tool('findPkgPath', package = 'CalibSvc') 
def exists(env):
    return 1;
