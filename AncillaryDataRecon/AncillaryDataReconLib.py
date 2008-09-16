# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['AncillaryDataRecon'])
    env.Tool('facilitiesLib')
    env.Tool('CalibDataLib')
    env.Tool('AncillaryDataEventLib')
    env.Tool('AncillaryDataUtilLib')
def exists(env):
    return 1;
