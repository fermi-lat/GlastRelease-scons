# $Header$
def generate(env, **kw):
     if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['AnalysisNtuple'])
        if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'AnalysisNtuple') 
     env.Tool('OverlayEventLib')
     env.Tool('EventLib')
     env.Tool('LdfEventLib')
     env.Tool('ntupleWriterSvcLib')
     env.Tool('OnboardFilterTdsLib')
     env.Tool('CalibDataLib')
     env.Tool('CalUtilLib')
     #env.Tool('RootIoLib')
     env.Tool('GlastSvcLib')
     env.Tool('AcdUtilLib')
     env.Tool('astroLib')
     env.Tool('facilitiesLib')
     env.Tool('addLibrary', library=env['gaudiLibs'])
     env.Tool('addLibrary', library=env['rootLibs'])
     env.Tool('addLibrary', library=env['clhepLibs'])
     if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
        env.Tool('findPkgPath', package = 'TkrUtil')
def exists(env):
    return 1;
