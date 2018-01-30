import subprocess
import numpy as np
import scipy.stats
import rebound


class CBSystem(object):
    '''Class to hold a circumbinary planet setup.'''

    def __init__(self,m1 = None, f1 = None, m2 = None, f2 = None,
                      ab = None, r1 = None, r2 = None, eb = None,
                      ib = np.pi/2., wb = None, fb = None,
                      ld11 = 0.85, ld21 = -0.09, ld12 = 0.85, ld22 = -0.09,
                      mp = 0.0 , ap = None, rp = None, ep = 0.0,
                      ip = np.pi/2., wp = None, fp = None, t0 = 0.0):
        '''Initialisation.
            
        Parameters
        ----------
        m1, m2, mp : float
            Mass of binary components and planet, in Solar masses.
        ab, eb, ib, fb, wb : float
            Orbital elements of binary, ab in au and radians for angles.
        ap, ep, ip, fp, wp : float
            Orbital elements of planet, ap in au and radians for angles.
        fl1, fl2 : float
            Fractional flux of binary components, planet is non-luminous.
        r1, r2, rp : float
            Radius of binary components and planet, in au.
        ld11, ld21, etc. : float
            Quadratic limb darkening coefficients, ld12 is first coeff
            for second body.
        t0 : float
            Time at which orbital elements apply, in days.
        '''
        self.m1 = m1
        self.f1 = f1
        self.m2 = m2
        self.f2 = f2
        self.ab = ab
        self.r1 = r1
        self.r2 = r2
        self.eb = eb
        self.ib = ib
        self.wb = wb
        self.fb = fb
        self.ld11 = ld11
        self.ld21 = ld21
        self.ld12 = ld12
        self.ld22 = ld22
        self.mp = mp
        self.ap = ap
        self.rp = rp
        self.ep = ep
        self.ip = ip
        self.wp = wp
        self.fp = fp
        self.t0 = t0


def heartbeat(sim):
    '''Heartbeat to compute transit times in python rebound.
    
    Observer is along positive z, so assume orbits are in x-z plane.
    '''
    
    global transittimes, lastt, lastdx1, lastdx2, lastdx12, lastdv1, lastdv2
    p = sim.contents.particles
    t = sim.contents.t
    
    # transit of primary
    thisdx1 = p[2].x -  p[0].x;
    thisdv1 = p[2].vx - p[0].vx;
    # condition is negative if different signs, want planet closer than star
    if (lastdx1*thisdx1 < 0) & (p[2].z > p[0].z):
        # interpolate transit time
        ttr = lastt  + (t-lastt)/(thisdx1-lastdx1)*(0.0-lastdx1)
        vtr = lastdv1+ (thisdv1-lastdv1)/(thisdx1-lastdx1)*(0.0-lastdx1)
        transittimes += ([20,ttr,np.abs(vtr)],)

    lastdx1 = thisdx1
    
    # transit of secondary
    thisdx2 = p[2].x -  p[1].x
    thisdv2 = p[2].vx - p[1].vx
    # condition is negative if different signs, want planet closer than star
    if (lastdx2*thisdx2 < 0) & (p[2].z > p[1].z):
        ttr = lastt  + (t-lastt)/(thisdx2-lastdx2)*(0.0-lastdx2);
        vtr = lastdv2+ (thisdv2-lastdv2)/(thisdx2-lastdx2)*(0.0-lastdx2);
        transittimes += ([21,ttr,np.abs(vtr)],)

    lastdx2 = thisdx2;

    # eclipses and occultations
    '''
    thisdx12 = thisdx2 - thisdx1   # distance between stars
    if lastdx12*thisdx12 < 0:
        if p[1].z > p[0].z:
            ttr = lastt+(t-lastt)/(thisdx12-lastdx12)*(0.0-lastdx12)
            transittimes +=  ([10,ttr,0.0],)
        else:
            ttr = lastt+(t-lastt)/(thisdx12-lastdx12)*(0.0-lastdx12)
            transittimes +=  ([1,ttr,0.0],)

    lastdx12 = thisdx12
    '''
    lastt = t


def reb_cb(cb,tmin=None,tmax=None):
    '''Return transit times for a circumbinary system.
    
    Parameters
    ----------
    cb : CBSys object
        Class holding system configuration.
    tmin, tmax : float
        Start and end times of simulation, in days.
    '''

    sim = rebound.Simulation()
    sim.t = 0.0
    sim.add(m = cb.m1)
    comp = rebound.Particle(simulation=sim, m=cb.m2, a=cb.ab, e=cb.eb,
                            inc=cb.ib, omega=cb.wb, f=cb.fb, Omega = 0.0)
    sim.add(comp)
    sim.move_to_com()
    planet = rebound.Particle(simulation=sim, m=cb.mp, a=cb.ap, e=cb.ep,
                              inc=cb.ip, omega=cb.wp, f=cb.fp, Omega=0.0)
    sim.add(planet)
    sim.move_to_com()

    # integrate to start time
    sim.integrate(tmax = (tmin - cb.t0)/365.25 * 2 * np.pi )

    # globals for heartbeat
    p = sim.particles
    global transittimes, lastt, lastdx1, lastdx2, lastdx12, lastdv1, lastdv2
    lastt = sim.t;
    lastdx1 = p[2].x -  p[0].x
    lastdx2 = p[2].x -  p[1].x
    lastdx12 = p[1].x -  p[0].x
    lastdv1 = p[2].vx - p[0].vx
    lastdv2 = p[2].vx - p[1].vx
    transittimes = ()

    # integrate to end
    sim.heartbeat = heartbeat
    sim.integrate(tmax = (tmax - cb.t0)/365.25 * 2 * np.pi )

    # get transit times and convert back to days
    tts = np.array(transittimes)
    tts[:,1] /= 2.0 * np.pi / 365.25
    tts[:,1] += cb.t0

    return tts


def reb_cb_c(cb, tmin=None, tmax=None,
             cb_path='/Users/grant/astro/code/github/rebound/examples/circumbinary'):
    '''As reb_cb, but using compiled c binary.'''

    # set times relative to tmin to avoid precision loss
    t0 = 0
    tmin_ = tmin - cb.t0
    tmax_ = tmax - cb.t0

    # set up and run the code
    cmd = [cb_path+'/rebound',
           '--t0='+str(t0),'--tmin='+str(tmin_),'--tmax='+str(tmax_),
           '--m1='+str(cb.m1),'--m2='+str(cb.m2),
           '--ab='+str(cb.ab),'--eb='+str(cb.eb),'--ib='+str(cb.ib),
           '--wb='+str(cb.wb),'--fb='+str(cb.fb),
           '--mp='+str(cb.mp),
           '--ap='+str(cb.ap),'--ep='+str(cb.ep),'--ip='+str(cb.ip),
           '--wp='+str(cb.wp),'--fp='+str(cb.fp)]

    x = subprocess.run(cmd,cwd=cb_path,stdout=subprocess.PIPE,check=True)

    # grab the output and make an array
    out = x.stdout.decode().rstrip().split('\n')
    tts = np.reshape( [s.split() for s in out], [len(out),3] ).astype(float)

    # convert back to days and add zero time
    tts[:,1] /= 2.0 * np.pi / 365.25
    tts[:,1] += cb.t0

    return tts


def convfm(f_in,ecc):
    '''Convert true anomaly to mean anomaly.
    
    Parameters
    ----------
    f_in : float
        True anomaly, in radians.
    ecc : float
        Eccentricity.
    '''
    tf2 = np.tan(0.5*f_in)
    fact = np.sqrt( (1.0+ecc)/(1.0-ecc) )
    bige = 2.0*np.arctan2(tf2,fact)
    bigm = bige - ecc*np.sin(bige)
    return bigm


def pd_cb(cb, times=None, filed=None, filet=None, cleanup=True,
          run='/Users/grant/Dropbox/astro/code/github/photodynam/photodynam'):
    '''Return flux from photodynam for a circumbinary system.
    
    Parameters
    ----------
    cb : CBSys object
        Class holding system configuration.
    times : list, tuple, or array
        Times to compute flux, in days.
    filed, filet : str
        File names for running photodynam.
    cleanup : bool, optional
        Remove files after running.
    run : str
        Path to photodynam.
    '''

    mcon = (365.25/2./np.pi)**2 # mass conversion factor, time is in days

    rstr = str(np.random.randint(1e5))
    if filed is None and filet is None:
        filed = '/tmp/pd1-'+rstr+'-.txt'
        filet = '/tmp/pd2-'+rstr+'-.txt'
    
    # write the dynamics file
    with open(filed,'w') as f:
        f.write('3 {}\n'.format(cb.t0))
        f.write('0.01 1e-16\n')                                 
        f.write('{} {} {}\n'.format(cb.m1/mcon, cb.m2/mcon, cb.mp/mcon))
        f.write('{} {} {}\n'.format(cb.r1, cb.r2, cb.rp))
        f.write('{} {} {}\n'.format(cb.f1, cb.f2, 0.0))
        f.write('{} {} {}\n'.format(cb.ld11, cb.ld12, 1.0))
        f.write('{} {} {}\n'.format(cb.ld21, cb.ld22, 0.0))
        # orbital elements, a, e, i, w, O=0.0, M
        f.write('{} {} {} {} 0.0 {}\n'.format(cb.ab, cb.eb, cb.ib, cb.wb,
                                              convfm( cb.fb, cb.eb )))
        f.write('{} {} {} {} 0.0 {}\n'.format(cb.ap, cb.ep, cb.ip, cb.wp,
                                              convfm( cb.fp, cb.ep )))
                 
    # write the times file, we want just flux as we know the times
    with open(filet,'w') as f:
        f.write('F\n')
        f.write(' '.join( map(str,times) ))
        
    # run the code and clean up
    x = subprocess.run([run,filed,filet],stdout=subprocess.PIPE,check=True)
    flux = np.genfromtxt(x.stdout.splitlines(),comments='#')

    if cleanup:
        subprocess.run(['rm',filed,filet])
    
    return flux


def stack(t, f, cb, window=1, event=20):
    '''Return a set of stacked light curve sections for a given system.
    
    Parameters
    ----------
    t : ndarray
        Array of times, in days.
    f : ndarray
        Array of fluxes or whatever is being stacked.
    window : float
        Width of region in days around transit times to include in stack.
    event : int
        Transit event to pick, 20 is planet-primary, 21 planet-secondary.
    '''

    all_tts = reb_cb_c(cb,tmin=np.min(t),tmax=np.max(t))

    # get the sections
    ok = all_tts[:,0] == event
    tts = all_tts[ok]
    stack = ()
    dates = ()
    for tt in all_tts[ok,1]:
        in_win = (-window < (t-tt)) & ((t-tt) < window)
        if np.any(in_win):
            stack += (f[in_win]/np.median(f[in_win]),)
            dates += (t[in_win] - tt,)

    return dates, stack


def stack_metric(ts, fs):
    '''Return a metric for the given stacked light curves.'''

    md,ed,n = scipy.stats.binned_statistic(np.hstack(ts), np.hstack(fs),
                                           'median', bins=50)

    return np.min(md)
