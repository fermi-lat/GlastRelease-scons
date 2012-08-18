#$Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'CalRecon') 

    #    env.Tool('addLibrary', library = ['CalRecon'])
    env.Tool('CalXtalResponseLib')
    env.Tool('TkrUtilLib')
    env.Tool('CalUtilLib')
    env.Tool('GuiSvcLib')
    env.Tool('GlastSvcLib')
    env.Tool('RootIoLib')
    #env.Tool('mcRootDataLib')
def exists(env):
    return 1;
