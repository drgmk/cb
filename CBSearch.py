
import numpy as np
import funcs
import pickle
import pylab as p
p.ion()

class Lightcurve():

    def __init__(self, time=None, flux=None, err=None, filepath=None, cut_outliers=False):

        if time:
            self.time = time
            self.flux = flux
            self.err = err
        elif filepath:
            lc = np.genfromtxt(filepath)
            self.time = lc[:,0]
            self.flux = lc[:,1]
            self.err = lc[:,2]
            if self.time[0] < 2000000:
                print('Adding Kepler Time offset of 2454833.')
                self.time += 2454833.
        if cut_outliers:
            self.cut_outliers()

    def cut_outliers(self,thresh=10,win=250):
        cutidx = funcs.outlier_cut(self.time,self.flux,thresh,win)
        self.time, self.flux, self.err = self.time[cutidx],self.flux[cutidx],self.err[cutidx]
    
    def cut_by_index(self, start, end):
        self.time = np.hstack((self.time[:start],self.time[end:]))
        self.flux = np.hstack((self.flux[:start],self.flux[end:]))
        self.err = np.hstack((self.err[:start],self.err[end:]))
        for win in self.windows:
            self.stattimes[win] =  np.hstack((self.stattimes[win][:start],self.stattimes[win][end:]))
            self.lcstat[win] =  np.hstack((self.lcstat[win][:start],self.lcstat[win][end:]))
            self.blurlcstat[win] =  np.hstack((self.blurlcstat[win][:start],self.blurlcstat[win][end:]))
            self.blurlcstat_norm[win] =  np.hstack((self.blurlcstat_norm[win][:start],self.blurlcstat_norm[win][end:]))
    #def cutphase
    
    def cut_transit(self,time,duration):
        start = np.searchsorted(self.time,time-duration/2.)
        end = np.searchsorted(self.time,time+duration/2.)
        self.time = np.hstack((self.time[:start],self.time[end:]))
        self.flux = np.hstack((self.flux[:start],self.flux[end:]))
        self.err = np.hstack((self.err[:start],self.err[end:]))
        #for win in self.windows:
        #    self.stattimes[win] =  np.hstack((self.stattimes[win][:start],self.stattimes[win][end:]))
        #    self.lcstat[win] =  np.hstack((self.lcstat[win][:start],self.lcstat[win][end:]))
        #    self.blurlcstat[win] =  np.hstack((self.blurlcstat[win][:start],self.blurlcstat[win][end:]))
        #    self.blurlcstat_norm[win] =  np.hstack((self.blurlcstat_norm[win][:start],self.blurlcstat_norm[win][end:]))



class CBSearch(Lightcurve):

    def __init__(self, lcfilepath, cut_outliers=False):
    
        Lightcurve.__init__(self, filepath=lcfilepath, cut_outliers=cut_outliers)
        self.tts_all = None
        self.periodogram = None
        self.per_archive = []
        self.windows = np.array([0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 
        				0.45, 0.5, 0.55, 0.6, 0.65, 0.7])
        self.minpoints = 3
        self.blurfactor = 1.
        self.minidx = None
        self.blurlcstat_norm = None
        
    def load_nbody(self,fpath1,fpath2=None):
    
        if self.periodogram:
            print('Resetting periodograms')
            self.periodogram = None
            self.per_archive = []
        
        with open(fpath1, 'rb') as handle:
            self.tts_all = {}
            self.tds_all = {}
            self.tts_all[0] = pickle.load(handle)
            self.tds_all[0] = pickle.load(handle)

        if fpath2:
            with open(fpath2, 'rb') as handle:
                self.tts_all[1] = pickle.load(handle)
                self.tds_all[1] = pickle.load(handle)

        self.set_keys()

    def set_keys(self):            
        self.ppset_keys = np.array(list(self.tts_all[0].keys()))
        self.ppset = self.ppset_keys.astype('float')
        self.fpset_keys = np.array(list(self.tts_all[0][list(self.tts_all[0].keys())[0]].keys()))
        self.fpset = self.fpset_keys.astype('float')
        
        if len(self.tts_all.keys())>1:
            self.ppset_2_keys = np.array(list(self.tts_all[1].keys()))
            self.ppset_2 = self.ppset_2_keys.astype('float')
            self.fpset_2_keys = np.array(list(self.tts_all[1][list(self.tts_all[1].keys())[0]].keys()))
            self.fpset_2 = self.fpset_2_keys.astype('float')

        #assert ((ppset==ppset_2).all()) 
               
    def drop_unstable(self):
        if self.tts_all:
            for pp in self.tts_all[0].keys():
                for fp in self.tts_all[0][pp].keys():
                    fail = False
                    if len(self.tts_all[0][pp][fp]) > np.ceil((self.time[-1]-self.time[0])/float(pp))+1:
                        fail = True
                    if len(self.tts_all.keys())>1:    
                        if len(self.tts_all[1][pp][fp]) > 2.5*np.ceil((self.time[-1]-self.time[0])/float(pp))+1:
                            fail = True
                    if fail:
                        self.tts_all[0][pp][fp] = []
                        if len(self.tts_all.keys())>1:
                            self.tts_all[1][pp][fp] = []
        else:
            print('No Nbody results loaded') 

    def drop_periods(self,ppcut=0.006):

        if self.periodogram:
            print('Resetting periodograms')
            self.periodogram = None
            self.per_archive = []
        
        if self.tts_all:
            pps_to_cut = []

            for pp in self.tts_all[0].keys():
                if 1./float(pp) > ppcut:
                    pps_to_cut.append(pp)
            for pp in pps_to_cut:
                self.tts_all[0].pop(pp)
                self.tds_all[0].pop(pp)
                if len(self.tts_all.keys())>1:
                    self.tts_all[1].pop(pp)
                    self.tds_all[1].pop(pp)
                
            self.set_keys() 
        else:
            print('No Nbody results loaded')    
                
    def make_lcstat(self, blurwindow=250):

        #window lightcurve with duration
        self.lcstat = {}
        self.blurlcstat = {}
        self.stattimes = {}
        self.blurantistat = {}
        
        print('Statisticking')
        for win in self.windows:
            self.lcstat[win], self.blurlcstat[win], self.stattimes[win], self.blurantistat[win] = funcs.running_mean_gaps(self.time,
            									     											self.flux-1,win,self.minpoints,self.blurfactor)
            self.blurlcstat[win] = self.blurlcstat[win] - np.median(self.blurlcstat[win])
            
        print('Normalising')
        self.blurlcstat_norm = funcs.normalise_stat(self.blurlcstat,self.lcstat,blurwindow)        #and norm
        self.blurantistat_norm = funcs.normalise_stat(self.blurantistat,self.lcstat,blurwindow)

        for win in self.windows:
            self.blurlcstat_norm[win][self.blurlcstat[win]>500] = 0.
            self.blurantistat_norm[win][self.blurlcstat[win]>500] = 0.
            self.blurantistat[win][self.blurlcstat[win]>500] = 0.
            self.stattimes[win][self.blurlcstat[win]>500] = 0.
            self.blurlcstat[win][self.blurlcstat[win]>500] = 0.
            self.lcstat[win][self.lcstat[win]>500] = 0.
            

    def make_periodogram(self,prisec=True,pertransit=False,target=0,statistic=None):
        if not statistic:
            if self.blurlcstat_norm:
                statistic = self.blurlcstat_norm
            else:
                print('No statistic provided and blurlcstat_norm not calculated')
                return 0
        
        if self.tts_all:
            if self.periodogram is not None:
                self.per_archive.append([int(self.prisec),int(self.pertransit),
            			    			 self.target,self.minidx,self.periodogram.copy()])

            self.prisec = prisec
            self.pertransit = pertransit
            self.target = target
        
            if prisec:
                if pertransit:
                    self.periodogram = funcs.make_periodogram_pertransit_prisec(
                    					self.tts_all[0],self.tts_all[1],self.tds_all[0],
                    					self.tds_all[1],self.time,self.ppset_keys,
                    					self.fpset_keys,self.windows,statistic)
                else:
                    self.periodogram = funcs.make_periodogram_prisec(
                    					self.tts_all[0],self.tts_all[1],self.tds_all[0],
                    					self.tds_all[1],self.time,self.ppset_keys,
                    					self.fpset_keys,self.windows,statistic)
            else:
                if pertransit:
                    self.periodogram = funcs.make_periodogram_pertransit(
                    					self.tts_all[target],self.tds_all[target],
                    					self.time,self.ppset_keys,self.fpset_keys,
                    					self.windows,statistic)
                else:
                    self.periodogram = funcs.make_periodogram(
                    					self.tts_all[target],self.tds_all[target],
                    					self.time,self.ppset_keys,self.fpset_keys,
                    					self.windows,statistic)
            self.minidx = np.unravel_index(np.argmin(self.periodogram),(len(self.ppset),len(self.fpset)))
        else:
            print('No Nbody results loaded')
            return 0

    def normalise_periodogram(self,statistic=None):
        '''
        Creates a periodogram using another statistic and subtracts it from the original.
        Not sure this is useful. Dividng would be intuitively more sensible, but has some issues.
        '''
        if not statistic:
            if self.blurlcstat_norm:
                statistic = self.blurlcstat_norm
            else:
                print('No statistic provided and blurlcstat_norm not calculated')
                return 0
                
        if self.periodogram is not None:
            if self.prisec:
                if self.pertransit:
                    self.antiperiodogram = funcs.make_periodogram_pertransit_prisec(
                    					self.tts_all[0],self.tts_all[1],self.tds_all[0],
                    					self.tds_all[1],self.time,self.ppset_keys,
                    					self.fpset_keys,self.windows,statistic)
                else:
                    self.antiperiodogram = funcs.make_periodogram_prisec(
                    					self.tts_all[0],self.tts_all[1],self.tds_all[0],
                    					self.tds_all[1],self.time,self.ppset_keys,
                    					self.fpset_keys,self.windows,statistic)
            else:
                if self.pertransit:
                    self.antiperiodogram = funcs.make_periodogram_pertransit(
                    					self.tts_all[self.target],self.tds_all[self.target],
                    					self.time,self.ppset_keys,self.fpset_keys,
                    					self.windows,statistic)
                else:
                    self.antiperiodogram = funcs.make_periodogram(
                    					self.tts_all[self.target],self.tds_all[self.target],
                    					self.time,self.ppset_keys,self.fpset_keys,
                    					self.windows,statistic)
            self.periodogram -= self.antiperiodogram
            self.minidx = np.unravel_index(np.argmin(self.periodogram),(len(self.ppset),len(self.fpset)))

        else:
            print('No periodogram calculated')
            return 0

        
        
    def set_from_archive(self,index):
    
        if self.periodogram:
            self.per_archive.append([int(self.prisec),int(self.pertransit),
            						 self.target,self.periodogram.copy()])

        if index < len(self.per_archive):
            self.prisec = self.per_archive[index][0]
            self.pertransit = self.per_archive[index][1]
            self.target = self.per_archive[index][2]
            self.minidx = self.per_archive[index][3]
            self.periodogram = self.per_archive[index][4]
        else:
            print('index too large, not enough archived periodograms')

    def plot_static(self, fppp0=None):
        if self.periodogram is None:
            print('No periodogram')
            return 0
        p.figure()
        palette = p.cm.viridis
        p.imshow(self.periodogram,origin='lower',
                extent=[np.min(self.fpset),np.max(self.fpset),1./np.min(self.ppset),1./np.max(self.ppset)],
                aspect='auto',cmap=palette)
        cbar = p.colorbar()
        cbar.set_label('Periodogram', rotation=270, labelpad=10)
        p.ylabel('Inverse Planet_period')
        p.xlabel('fp')
        if fppp0 is not None:
            p.plot(fppp0[0],1./fppp0[1],'rx')
        p.plot(self.fpset[self.minidx[1]],1./self.ppset[self.minidx[0]],'bx')
        print(np.min(self.periodogram))
        print((np.min(self.periodogram)-np.mean(self.periodogram))/np.std(self.periodogram))
        

    def plot_diag(self, fppp0=None):
        if self.periodogram is None:
            print('No periodogram')
            return 0
        fig = p.figure(1)
        p.clf()
        palette = p.cm.viridis
        p.imshow(self.periodogram,origin='lower',
                    extent=[np.min(self.fpset),np.max(self.fpset),1./np.min(self.ppset),1./np.max(self.ppset)],
                    aspect='auto',cmap=palette)
        cbar = p.colorbar()
        cbar.set_label('Periodogram', rotation=270, labelpad=10)
        p.ylabel('Inverse Planet_period')
        p.xlabel('fp')
        if fppp0 is not None:
            p.plot(fppp0[0],1./fppp0[1],'rx')
        p.plot(self.fpset[self.minidx[1]],1./self.ppset[self.minidx[0]],'bx')

        #global ix, iy
        ix, iy = [0], [1000]
        
        def onclick(event):
            p.figure(1)
            tb = p.get_current_fig_manager().toolbar 
            if event.button==1 and event.inaxes and tb.mode == '': #check not zooming etc
                #global ix, iy
                cut = False
                if np.abs(ix[0]-event.xdata)<0.01 and np.abs(iy[0]-event.ydata)<0.01:
                    cut = True

                ix[0], iy[0] = event.xdata, event.ydata
                #tx = 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % (event.button, event.x, event.y, event.xdata, event.ydata)
                #text.set_text(tx)
                print(ix[0],iy[0])
            
                #identify window - find nearest tts, identify one with the maximum effect
                self.minidx = [np.argmin(np.abs(self.ppset-1./iy[0])),np.argmin(np.abs(self.fpset-ix[0]))]
                nearpp = self.ppset_keys[self.minidx[0]]
                nearfp = self.fpset_keys[self.minidx[1]]
             
                tts = self.tts_all[0][nearpp][nearfp]
                tds = self.tds_all[0][nearpp][nearfp]
                if len(self.tts_all.keys())>1:
                    tts = np.hstack((tts,self.tts_all[1][nearpp][nearfp]))
                    tds = np.hstack((tds,self.tds_all[1][nearpp][nearfp]))
            
                maxtime = 0
                maxstat = 0
                toplot = []
                for t in range(len(tts)):
                    window = self.windows[np.argmin(np.abs(self.windows-tds[t]))]
                    timeindex = np.searchsorted(self.time,tts[t])
                    if timeindex < len(self.blurlcstat_norm[window]):
                        if self.time[timeindex] - tts[t] < window:
                            toplot.append([timeindex,window])
                            if self.blurlcstat_norm[window][timeindex] < maxstat:
                                maxstat = self.blurlcstat_norm[window][timeindex]
                                maxtime = self.stattimes[window][timeindex]
                                maxdur = window     
                
                if cut:
                    start = np.searchsorted(self.time,maxtime-maxdur/2.)
                    end = np.searchsorted(self.time,maxtime+maxdur/2.)
                    self.cut_by_index(start,end)
                    print('Section Removed')
                else:
                    print('Pixel set. Click near pixel again to remove this section.')

                p.figure(2)
                p.clf()
                f, ax = p.subplots(len(toplot),num=2,figsize=(6.4,2.*len(toplot)))
                    
                for t,transit in enumerate(toplot):
                    window = transit[1]
                    ttime = self.stattimes[window][transit[0]]
                    start = np.searchsorted(self.time,ttime-2.)
                    end = np.searchsorted(self.time,ttime+2.)
                    fluxlevel = np.mean(self.flux[start:end])-3*np.std(self.flux[start:end])
                    if ttime == maxtime:
                        ax[t].plot(self.time[start:end],self.flux[start:end],'r.')
                    else:
                        ax[t].plot(self.time[start:end],self.flux[start:end],'b.')
                    ax[t].plot([ttime-window/2.,ttime+window/2.],[fluxlevel,fluxlevel],'m.-')
                    
                p.pause(0.1)

                p.figure(1) #set current fig to the main periodogram, for toolbar check
                    
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
                
    def plot_transits(self,win=None):
        if not win:
            win = self.windows[4]
        if self.minidx:
            p.figure()
            p.plot(self.time,self.flux,'b.')
            p.plot(self.time,self.lcstat[win],'r.-')
            p.plot(self.time,self.blurlcstat[win],'g.-')

            pp = self.ppset_keys[self.minidx[0]]
            fp = self.fpset_keys[self.minidx[1]]
            #pp = ppset_keys[29]
            #fp = fpset_keys[572]

            tts = self.tts_all[1][pp][fp]
            tds = self.tds_all[1][pp][fp]
            if len(self.tts_all.keys())>1:
                tts = np.hstack((tts,self.tts_all[1][pp][fp]))
                tds = np.hstack((tds,self.tds_all[1][pp][fp]))

            for transit,dur in zip(tts,tds):
                p.plot(transit,0.9985,'m.')
                #p.plot([transit-dur/2.,transit+dur/2.],[0.999,0.999])

            p.ylim(0.998,1.002)    
        
    def plot_pfold(self,ndurs=11.,statistic=None):
        if not statistic:
            if self.blurlcstat_norm:
                statistic = self.blurlcstat_norm
            else:
                print('No statistic provided and blurlcstat_norm not calculated')
                return 0
                
        if self.minidx:
            pp = self.ppset_keys[self.minidx[0]]
            fp = self.fpset_keys[self.minidx[1]]
            #pp = ppset_keys[29]
            #fp = fpset_keys[572]

            tts = self.tts_all[0][pp][fp]
            tds = self.tds_all[0][pp][fp]
            if len(self.tts_all.keys())>1:
                tts = np.hstack((tts,self.tts_all[1][pp][fp]))
                tds = np.hstack((tds,self.tds_all[1][pp][fp]))

            p.figure()
            #intransit = np.zeros(0)

            for transit,dur in zip(tts,tds):
                time_window, flux_window, timescale = funcs.extract_transit_window_sourcetimes(
                										transit,dur,self.time,self.flux,
                										self.windows,self.stattimes,statistic,float(ndurs)) 
                start = np.searchsorted(timescale,-1/float(ndurs))
                end = np.searchsorted(timescale,1/float(ndurs))
                #intransit = np.hstack((intransit,flux_window[start:end]))

                if len(time_window)>0:
                    p.plot(timescale,flux_window,'m.')
                    p.plot([-1/float(ndurs),-1/float(ndurs)],[np.min(flux_window),np.max(flux_window)],'r--')
                    p.plot([1/float(ndurs),1/float(ndurs)],[np.min(flux_window),np.max(flux_window)],'r--')
    
    def plot_hist(self,bins=100):
        if self.minidx:
            p.figure()
            p.hist(self.periodogram.flatten(),bins=bins)
            p.yscale('log', nonposy='clip')    


