'''Run n-body simulation of planet around binary to get timings.'''

import pickle
import copy
import numpy as np
import matplotlib.pyplot as plt
import itertools
from multiprocessing import Pool

import funcs
import known_system_params

#def run_cb(p):
def run_cb(pp, fp, cb, transited_comp, t_min, t_max, timing_precision):
    '''Wrapper to run n-body.'''

    cb.ap = ( (cb.m1+cb.m2) * (pp/365.25)**2 )**(1./3.)
    cb.fp = fp    
    try:
        return funcs.reb_cb_dvm(cb, transited_comp, 2, tmin=t_min,
                                tmax=t_max,
                                timing_precision=timing_precision)
    except:
        return np.array([]),np.array([])


def run_nbody_grid(system, n_period_bins=None, n_f_bins=700,
                   transited_comp=0, t_min=2454963., t_max=2456424.,
                   min_a_ab=3., degrade=True, lc_file=None,
                   n_processes=8, out_dir='../nbody_data/'):

    # get a system
    if system == 'k16':
        cb = known_system_params.k16
    elif system == 'k34':
        cb = known_system_params.k34
    elif system == 'k35':
        cb = known_system_params.k35
    elif system == 'k47b':
        cb = known_system_params.k47b
    elif system == 'k47c':
        cb = known_system_params.k47c
    elif system == 'k64':
        cb = known_system_params.k64
    elif system == 'k413':
        cb = known_system_params.k413
    elif system == 'k453':
        cb = known_system_params.k453
    elif system == 'k1647':
        cb = known_system_params.k1647

    print('system: {}'.format(system))

    # we will modify this later, so ensure we have a copy of the original
    cb = copy.deepcopy(cb)

    # binary period, in days
    p_bin = (cb.ab**3/(cb.m1+cb.m2))**(1./2.)*365.25

    # time precision of simulation
    timing_precision = 30./86400./365.25 * (2*np.pi)  # in years/2pi

    # load lightcurve and check stack looks OK
    if lc_file is not None:
        lc = np.genfromtxt(lc_file)
        time = lc[:,0]+2454833
        flux = lc[:,1]
        err = lc[:,2]

        # and check it stacks right
        window = 1.0
            
        fig,ax = plt.subplots(2, figsize=(9.5,5), sharex=True)

        times, stack = funcs.stack(time, flux, cb)
        for t1,s1 in zip(times,stack):
            ax[0].plot(t1, s1)

        times, stack = funcs.stack(time, flux, cb, event=21)
        for t1,s1 in zip(times,stack):
            ax[1].plot(t1, s1)
            
        ax[0].set_title('transit of primary for nominal fit')
        ax[1].set_title('transit of secondary for nominal fit')

        fig.savefig(out_dir+system+'-'+str(transited_comp)+'-nbody-stack.png')

    # degrading our knowledge (assume we know eb and wb):
    # radii left as we know them - this will change transit durations,
    # and whether some marginal transits happen.
    if degrade:
        cb.ib = np.pi/2.
        cb.ip = np.pi/2.
        cb.Wp = 0.0
        cb.ep = 0.0
        #cb.wp left as is because it defines f_p
        cb.mp = 0.0

    #define range over which to scan
    per_s_limit = p_bin*min_a_ab**1.5
    if n_period_bins is None:
        n_period_bins = int( (1/per_s_limit - 1/500.) / (2.5e-6) )

    invppset = np.linspace(1/500.,1/per_s_limit,n_period_bins)
    ppset = 1./invppset[::-1] #uniform in frequency
    ppname = [str(pp)[:6] for pp in ppset]

    fpset = np.linspace(0,2*np.pi,n_f_bins)
    fpname = [str(fp)[:6] for fp in fpset]

    print('n period bins: {}'.format(len(ppset)))
    print('n anomaly bins: {}'.format(len(fpset)))

    # loop over periods until no stable orbits
    tts_all = {}
    tds_all = {}
    for ps,pns in zip(ppset, ppname):
        
        # cartesian product -> all params and names in a big list
        par = itertools.product([ps],fpset)
        par_all = [p+(copy.deepcopy(cb), transited_comp,
                      t_min, t_max, timing_precision) for p in par]
        parname = itertools.product([pns],fpname)

        # run it
        pool = Pool(processes=n_processes)
        out = pool.starmap(run_cb, par_all)
        pool.close()

        tts_all[pns] = {}
        tds_all[pns] = {}
       
        for l,n in zip(out,parname):
            tt, td = l
            pn, fn = n
            tts_all[pn][fn] = tt[td>0]
            tds_all[pn][fn] = td[td>0]

    with open(out_dir+system+'-'+str(transited_comp)+'-nbody.pkl','wb') as f:
        pickle.dump(tts_all,f)
        pickle.dump(tds_all,f)


# run them
comp = 0
#run_nbody_grid('k16',lc_file='Lightcurves/K16_clean_noecl.txt',
#               transited_comp=comp)
#run_nbody_grid('k34',lc_file='Lightcurves/K34_clean_noecl.txt',
#               transited_comp=comp)
#run_nbody_grid('k35',lc_file='Lightcurves/K35_clean_noecl.txt',
#               transited_comp=comp)
run_nbody_grid('k47b',lc_file='Lightcurves/K47_clean_noecl.txt',
               transited_comp=comp)
