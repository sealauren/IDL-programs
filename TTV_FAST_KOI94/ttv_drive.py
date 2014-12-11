G = 0.000295994511 # G in units of AU, day

def mask_rm_transit(rvdata, t_rm):
    rvdata2 = rvdata.copy()
    rvdata2 = rvdata2[(min(t_rm) > rvdata2.time ) | (rvdata2.time > max(t_rm))]
    return rvdata2
    
def grow_text(bigbodies):
	from numpy import arange
	npl = len(bigbodies.Period)
	text=str(G)+'\n'
	text+=str(bigbodies.Mstar[0])+"\n"
	for j in arange(npl):
		text+=str(bigbodies.Mplanet[j])+'\n'
		text+=str(bigbodies.Period[j])+' '
		text+=str(bigbodies.E[j])+' '
		text+=str(bigbodies.I[j])+' '
		text+=str(bigbodies.LongNode[j])+' '
		text+=str(bigbodies.Argument[j])+' '
		text+=str(bigbodies.MeanAnomaly[j])+' '+"\n"
	return text
    
def write_input_file(text, path="/Users/lweiss/TTV_FAST_KOI94/", fileid='1'):
    text_file = open(path+"test.in."+fileid, "w")
    text_file.write(text)
    print text
    text_file.close()
    
def write_setup_file(npl,T0,path="/Users/lweiss/TTV_FAST_KOI94/", fileid='1'):
    setup_file = open(path+"setup_file", "w")
    setup_text = 'test.in.'+str(fileid)+'\n'+str(T0)+'\n'+'0.10\n'+str(float(T0)+2000)+'\n'+str(npl)+'\n'+'0\n'
    setup_file.write(setup_text)
    print setup_text
    setup_file.close()
    
def write_rv_file2(T0, path="/Users/lweiss/TTV_FAST_KOI94/"):
    from numpy import arange
    RV_times = arange(T0,T0+2000.,0.1)
    file = open(path+"RV_file2", "w")
    for t in RV_times:
        file.write('     '+str(t)+'\n')
    file.close()
    
def run_TTV_FAST(fileid):
    from subprocess import Popen, PIPE, STDOUT
    cmd1 = "./run_TTVFast setup_file Times RV_file RV_out"
    cmd2 = "./run_TTVFast setup_file Times RV_file2 RV_out2"
    final = Popen("{}; {}".format(cmd1, cmd2), shell=True, stdin=PIPE, 
              stdout=PIPE, stderr=STDOUT, close_fds=True)
    stdout, nothing = final.communicate()
    log = open('/log', 'w')
    log.write(stdout)
    log.close()
    return
    
def random_draw_MplPEIOwM(big_bodies):
    Mpl = np.random.normal(scale=1.0,size=4)*big_bodies.Mplanet*0.01 + big_bodies.Mplanet
    P = np.random.normal(scale=1.0,size=4)*big_bodies.Period*0.01 + big_bodies.Period
    E = np.random.normal(scale=1.0,size=4)*big_bodies.E*0.1 + big_bodies.E
    I = np.random.normal(scale=1.0,size=4)*big_bodies.I*0.01 + big_bodies.I
    O = (np.random.normal(scale=1.0,size=4)*big_bodies.LongNode*0.1 + big_bodies.LongNode) % 360.
    w = (np.random.normal(scale=1.0,size=4)*0.1 + big_bodies.Argument) % 360.
    M = (np.random.normal(scale=1.0,size=4)*big_bodies.MeanAnomaly*0.01 + big_bodies.MeanAnomaly) % 360.
    return(Mpl,P,E,I,O,w,M)
    
def read_results():
    Times = pd.read_csv('Times', header=None, names=['Planet','Epoch','Time','RSKY_AU','VSKY_AU'],sep=' ',usecols=[0,1,2,3,4])
    RV_out = pd.read_csv('RV_out', header=None, names=['time', 'rv'],sep=' ',usecols=[0,1])
    RV_out2 = pd.read_csv('RV_out2', header=None, names=['time', 'rv'],sep=' ',usecols=[0,1])
    offset = 0.0
    AU = 149597870700.0 #in meters
    day = 3600.0*24.0 #in seconds
    RV_out.rv = RV_out.rv*AU/day+offset
    RV_out2.rv = RV_out2.rv*AU/day+offset
    return Times, RV_out, RV_out2
    
def rv_resid(rvdata, RV_out):
    return np.array(RV_out.rv) - np.array(rvdata.rv)
    
def tt_resid(ttvdata, Times, Planet=2, verbose=False):
    dtarray = []
    for i in arange(len(ttvdata.time)):
        if verbose==True: print "TT (data):", ttvdata.time[i], type(ttvdata.time[i])
        deltat = min(abs(Times[Times.Planet==Planet].Time - ttvdata.time[i]))
        if verbose==True: print "Delta T:", deltat
        dtarray.append(deltat)
    dtarray = np.array(dtarray)
    return dtarray
    
def write_pars(npl,P0,tp0=None, ecc0=None, om0=None, K0=None, tt0=None):
    import numpy as np
    from lmfit import minimize, Parameters, Parameter, report_fit
    assert len(P0) == npl
    print "Writing parameters for %s planets..." %npl
    if om0==None: om0 = np.zeros(npl) + 90.
    if ecc0==None: ecc0 = np.zeros(npl) + 0.1
    if tp0==None: tp0 = np.zeros(npl) + min(tarr)
    if K0==None: K0 = np.zeros(npl) + 5.0
    pars = Parameters()
    pdict = {0:'b', 1:'c', 2:'d', 3:'e', 4:'f', 5:'g'}
    for i in arange(npl):
        pars.add_many(('P_'+str(pdict[i])    , P0[i]  , False,  None,  None,  None),
                    ('tt_'+str(pdict[i]), tt0[i], False, None, None, None),
                    ('om_'+str(pdict[i]), om0[i]       , True, 0., 360.,  None),
                    ('logK_'+str(pdict[i])  , np.log(K0[i])  , True, None, None,  None),
                    ('tp_'+str(pdict[i])    , tp0[i]  , True,  None,  None,  
                     'compute_tp(eval("tt_"+pdict[0]),0.,eval("P_"+pdict[0]),eval("om_"+pdict[0]))'),
                    ('ecc_'+str(pdict[i]), ecc0[i]     , True,  0., .5,  None))
    pars.add_many(('gamma' , 0.0 , True, None, None, None),
                   ('dvdt'  , 0.    , False, None, None, None),
                   ('curve', 0., False, None, None, None),
                   ('tstar',tp0[0], False, None, None))
    return pars