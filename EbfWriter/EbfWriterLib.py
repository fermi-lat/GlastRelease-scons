#$Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['EbfWriter'])
    env.Tool('CalXtalResponseLib')
    env.Tool('TriggerLib')
    env.Tool('CalUtilLib')
    env.Tool('FluxSvcLib')
    env.Tool('fluxLib')
    env.Tool('EventLib')
    env.Tool('GlastSvcLib')
    env.Tool('astroLib')
def exists(env):
    return 1;
