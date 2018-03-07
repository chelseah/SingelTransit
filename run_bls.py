#/usr/bin/env python
import matplotlib
import numpy as np
import scipy as sp
from dataio import readcolumn
from matplotlib import pyplot as plt

def phasefold(time, period, epoch):
    # calculate the phase given a period of a transit, the center of transit is phase 0
    ftime = np.zeros(len(time))
    ftime = (time-epoch-0.5*period)/period-((time-epoch-0.5*period)/period).astype(int)
    ind=ftime<0
    ftime[ind]+=1
    return ftime

def cal_bls(phase, mag, q): 
    # compute the BLS statistic given phase, magnitude and q (fraction of transit given period), see kovacs (2002?) for all the definetions
    w = 1./len(phase)
    mavg = np.mean(mag)
    w2sum = w
    m2avg = np.std(mag)**2.
    intran = (phase>0.5-q/2.)*(phase<0.5+q/2.)
    r = len(mag[intran])
    s = sum(mag[intran]-mavg)*w
    r2 = len(mag[intran])*w**2.
    s2 = sum((mag[intran]-mavg)**2.)*w
    if r == 0:
        return [np.nan, np.nan, np.nan]
    L = s/r+mavg
    H = s/(r-1.0)+mavg
    dip = L-H
    Lvar = (r*s2-s**2.)/(r**2-r2)
    Hvar = ((m2avg-s2)*(1.-r)-s**2.)/((1-r)**2.-(w2sum-r2))
    #print Lvar*r2/r**2+Hvar*(w2sum-r2), (1.-r)**2. 
    dipsig = np.sqrt(np.abs(Lvar*r2/r**2+Hvar*(w2sum-r2)/(1.-r)**2.))
    dsp = abs(dip/dipsig)
    return [dsp,dip,dipsig]

def bls_fixed_epoch_duration(time, mag, epoch, Tdur):
    # compute the BLS spectra of a known epoch (searching for period) 
    
     
    periodgrid = 1./np.linspace(1./30., 1./0.5, 26000.)
    qvararr = np.zeros(len(periodgrid))
    dsparr = np.zeros(len(periodgrid))
    diparr = np.zeros(len(periodgrid))
    dipsigarr = np.zeros(len(periodgrid))
    for i in xrange(len(periodgrid)):
        phase = phasefold(time, periodgrid[i], epoch)
        #print epoch, max(time), min(time), max(phase), min(phase)
        #break 
        qvar = Tdur/periodgrid[i]
        dsp, dip, dipsig = cal_bls(phase, mag, qvar)
        
        # print i, periodgrid[i], qvar, dsp, dip, dipsig
        qvararr[i] = qvar
        dsparr[i] = dsp
        diparr[i] = dip
        dipsigarr[i] = dipsig
    return [periodgrid, qvararr, dsparr, diparr, dipsigarr]
def bls_fixed_period_duration(time, mag, period, Tdur):
    # compute the BLS spectra of a known epoch (searching for period) 
    
     
    epochgrid = np.linspace(0, 1, 10000)
    dsparr = np.zeros(len(epochgrid))
    diparr = np.zeros(len(epochgrid))
    dipsigarr = np.zeros(len(epochgrid))

    for i in xrange(len(epochgrid)):
        phase = phasefold(time, period, time[0]+period*epochgrid[i])
        #print epoch, max(time), min(time), max(phase), min(phase)
        #break 
        qvar = Tdur/period
        dsp, dip, dipsig = cal_bls(phase, mag, qvar)
        dsparr[i] = dsp
        diparr[i] = dip
        dipsigarr[i] = dipsig
    return [epochgrid, dsparr, diparr, dipsigarr]




def example_period_search():
    infile = "Kepler-6b_ltf.lc" 
    time = []; readcolumn(time, 1, infile); time = np.array(time)
    flux = []; readcolumn(flux, 3, infile); flux = np.array(flux)
    epoch = 954.48636
    Tdur = 1./24. 
    periodgrid, qvararry, dsparr, diparr, dipsigarr = bls_fixed_epoch_duration(time, flux, epoch, Tdur)

    plt.plot(periodgrid, dsparr, 'k')
    plt.show()



def example_epoch_search():
    #infile = "Kepler-6b_ltf.lc" 
    infile = "/Volumes/Portabella/Data/TESS/Stamp/ASCIILC/270826460.rlc" 
    time = []; readcolumn(time, 2, infile); time = np.array(time)
    flux = []; readcolumn(flux, 19, infile); flux = np.array(flux)
    time/=48.
    import scipy as sp
    from scipy.signal import medfilt
    rawflux = flux
    base = medfilt(flux, 1055)
    #plt.plot(time, base, '.')
    #plt.show()
    flux = flux/base


    period = time[-1]-time[0] 
    Tdur = 1./24. 
    epochgrid, dsparr, diparr, dipsigarr = bls_fixed_period_duration(time, flux, period, Tdur)

    sigma = np.nanmedian(np.abs(dsparr-np.nanmedian(dsparr)))
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.plot(epochgrid*period/15., dsparr/sigma, 'k')
    ax1.set_xlabel("Time (day)")
    ax1.set_ylabel("BLS SNR")
    #plt.plot(epochgrid, dipsigarr, 'k')
    
    ax2 = fig.add_subplot(212)
    ax2.plot((time-time[0])/15., np.mean(rawflux)-rawflux, '.')
    ax2.set_xlabel("Time (day)")
    ax2.set_ylabel("Relative Flux")
    plt.show()
    


if __name__=='__main__':
   #example_period_search() 
   example_epoch_search()
