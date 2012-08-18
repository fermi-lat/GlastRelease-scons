#$Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['geomrep'])
        if env['PLATFORM'] == "win32" and env.get('CONTAINERNAME','') == 'GlastRelease':
            env.Tool('findPkgPath', package = 'geomrep') 
            env.Tool('findPkgPath', package = 'geometry') 
            env.Tool('findPkgPath', package = 'gui') 
    env.Tool('geometryLib')
    env.Tool('guiLib')
    if kw.get('incsOnly', 0) == 1: 
        env.Tool('findPkgPath', package = 'gui')
        env.Tool('findPkgPath', package = 'geometry')

def exists(env):
    return 1;
