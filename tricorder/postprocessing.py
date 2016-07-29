from astropy.io import fits
import treecorr
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json
from scipy.stats import binned_statistic

# --------------------------
#Eventually should be object-oriented following this skeleton:

class NNNPlotter (object):
    def __init__(self,runname):
        config_fname = outdir+runname+'.config'
        try: 
            with open(config_fname) as f:
                configdict = json.loads(f.read())
        except IOError:
            print 'results for '+runname+' not found in '+outdir

    def load(self):
        [ddd,ddr,drd,rdd,drr,rdr,rrd,rrr] = 8*[treecorr.NNNCorrelation(config=self.config)]
        for nnn in [ddd,ddr,drd,rdd,drr,rdr,rrd,rrr]:
            pass

# -----------------------------
        
def computeAngularBins(r,u,v,collapsed=False):
    d2 = r
    d3 = u*r
    d1 = v*d3+d2
    #law of cosines
    if not collapsed:
        cosine = (d2**2 + d3**2 - d1**2)/(2*d2*d3+1e-9)
    else: 
        cosine = (d1**2 + d3**2 - d2**2)/(2*d1*d3+1e-9)
    bins = np.arccos(cosine)/np.pi*180
    return bins

def computeXvsAngle(var):
    nbins = 8
    out1,bins1,_ = binned_statistic(elongated_angles[np.where(isRightSize & isElongated)],var[np.where(isRightSize & isElongated)],bins=nbins,statistic='mean')
    out2,bins2,_ = binned_statistic(collapsed_angles[np.where(isRightSize & isCollapsed)],var[np.where(isRightSize & isCollapsed)],bins=nbins,statistic='mean')
    full_var = np.concatenate((out2,out1))
    #make edges centers
    bins1 += (bins1[1]-bins1[0])/2
    bins2 += (bins2[1]-bins2[0])/2
    full_bins = np.concatenate((bins2[:-1],bins1[:-1]))
    return full_var, full_bins

output_dir = '/nfs/slac/g/ki/ki19/des/mbaumer/3pt_runs/from_sherlock/'
runname = 'y1_allz'

#load data
with open(output_dir+runname+'.config') as f:
    config = json.loads(f.read())
ddd = treecorr.NNNCorrelation(config=config)
drr = treecorr.NNNCorrelation(config=config)
rdr = treecorr.NNNCorrelation(config=config)
rrd = treecorr.NNNCorrelation(config=config)
ddr = treecorr.NNNCorrelation(config=config)
drd = treecorr.NNNCorrelation(config=config)
rdd = treecorr.NNNCorrelation(config=config)
rrr = treecorr.NNNCorrelation(config=config)
ddd.read(output_dir+runname+'ddd.out',file_type='ASCII')
drr.read(output_dir+runname+'drr.out',file_type='ASCII')
rdr.read(output_dir+runname+'rdr.out',file_type='ASCII')
rrd.read(output_dir+runname+'rrd.out',file_type='ASCII')
ddr.read(output_dir+runname+'ddr.out',file_type='ASCII')
drd.read(output_dir+runname+'drd.out',file_type='ASCII')
rdd.read(output_dir+runname+'rdd.out',file_type='ASCII')
rrr.read(output_dir+runname+'rrr.out',file_type='ASCII')

data = fits.getdata(output_dir+'redmagic_Y1_sims_data.fits')
randoms = fits.getdata(output_dir+'redmagic_Y1_sims_5x_randoms.fits')

data = data[((data['ZREDMAGIC'] > config['min_z']) & (data['ZREDMAGIC'] < config['max_z']))]
data = data[((data['RA'] > config['min_ra']) & (data['RA'] < config['max_ra']))]
data = data[((data['DEC'] > config['min_dec']) & (data['DEC'] < config['max_dec']))]

randoms = randoms[((randoms['Z'] > config['min_z']) & (randoms['Z'] < config['max_z']))]
randoms = randoms[((randoms['RA'] > config['min_ra']) & (randoms['RA'] < config['max_ra']))]
randoms = randoms[((randoms['DEC'] > config['min_dec']) & (randoms['DEC'] < config['max_dec']))]

cat = treecorr.Catalog(ra=data['RA'], dec=data['DEC'], ra_units='degrees', dec_units='degrees')
random_cat = treecorr.Catalog(ra=randoms['RA'], dec=randoms['DEC'], ra_units='degrees', dec_units='degrees')

#treecorr hack
ddd.tot = float(cat.nobj)**3/6
drr.tot = float(cat.nobj)*float(random_cat.nobj)**2/6
rdr.tot = float(cat.nobj)*float(random_cat.nobj)**2/6
rrd.tot = float(cat.nobj)*float(random_cat.nobj)**2/6
ddr.tot = float(cat.nobj)**2*float(random_cat.nobj)/6
drd.tot = float(cat.nobj)**2*float(random_cat.nobj)/6
rdd.tot = float(cat.nobj)**2*float(random_cat.nobj)/6
rrr.tot = float(random_cat.nobj)**3/6

zeta, varzeta = ddd.calculateZeta(drr=drr,rdr=rdr,rrd=rrd,ddr=ddr,drd=drd,rdd=rdd,rrr=rrr)

collapsed_angles = computeAngularBins(np.exp(ddd.logr),ddd.u,ddd.v,collapsed=True)
elongated_angles = computeAngularBins(np.exp(ddd.logr),ddd.u,ddd.v,collapsed=False)

isRightSize = ((ddd.meand3 > 7.8) & (ddd.meand3 < 8.2))#((ddd.meand1 > 5) & (ddd.meand1 < 19) & (ddd.meand2 > 5) & (ddd.meand2 < 19) & (ddd.meand3 > 5) & (ddd.meand3 < 7))
#isCollapsed = ((ddd.meand1 > 10) & (ddd.meand1 < 14))
#isElongated = ((ddd.meand2 > 11) & (ddd.meand2 < 13))
isCollapsed = ((ddd.meand1/ddd.meand3 > 1.75) & (ddd.meand1/ddd.meand3 < 2.25) & (ddd.v > 0))
isElongated = ((ddd.u > .45) & (ddd.u < .55) & (ddd.v > 0))
isEquilateral = (((np.abs(ddd.meand1 - ddd.meand2) < 1) & (np.abs(ddd.meand2 - ddd.meand3) < 1) & (np.abs(ddd.meand1 - ddd.meand3) < 1)))
isInSample = (isRightSize & (isCollapsed | isElongated))
sampled1 = ddd.meand1[np.where(isInSample)]
sampled2 = ddd.meand2[np.where(isInSample)]
sampled3 = ddd.meand3[np.where(isInSample)]
collapsedV = ddd.v[np.where(isRightSize & isCollapsed)]
collapsedU = ddd.u[np.where(isRightSize & isCollapsed)]
elongatedV = ddd.v[np.where(isRightSize & isElongated)]
elongatedU = ddd.u[np.where(isRightSize & isElongated)]
equilateralR = np.exp(ddd.logr[np.where(isEquilateral)])
equilateralNtri = ddd.ntri[np.where(isEquilateral)]

zetabins,bins = computeXvsAngle(zeta,)
d1bins, bins = computeXvsAngle(ddd.meand1)
d2bins, bins = computeXvsAngle(ddd.meand2)
d3bins, bins = computeXvsAngle(ddd.meand3)

fig=plt.figure()
plt.plot(bins,zetabins,'.')
plt.xlabel(r'$\theta$ (degrees)')
plt.ylabel('Zeta')
plt.title('Zeta vs. angle')
fig.savefig(runname+'_ZetavsAngle.png')

fig=plt.figure()
plt.plot(bins,d1bins,'.')
plt.xlabel(r'$\theta$ (degrees)')
plt.ylabel('d1')
plt.title('d1 vs. angle')
fig.savefig(runname+'_d1vsAngle.png')

fig=plt.figure()
plt.plot(bins,d2bins,'.')
plt.xlabel(r'$\theta$ (degrees)')
plt.ylabel('d2')
plt.title('d2 vs. angle')
fig.savefig(runname+'_d2vsAngle.png')

fig=plt.figure()
plt.plot(bins,d3bins,'.')
plt.xlabel(r'$\theta$ (degrees)')
plt.ylabel('d3')
plt.title('d3 vs. angle')
fig.savefig(runname+'_d3vsAngle.png')

dd = treecorr.NNCorrelation(config=config,nbins=10)
rr = treecorr.NNCorrelation(config=config,nbins=10)
dd.process(cat)
rr.process(random_cat)

xi, varxi = dd.calculateXi(rr=rr)
plt.loglog(np.exp(dd.logr),xi)
coeffs = np.polyfit(dd.logr,np.log(xi),deg=1)
poly = np.poly1d(coeffs)
yfit = lambda x: np.exp(poly(np.log(x)))

xi1 = yfit(d1bins)
xi2 = yfit(d2bins)
xi3 = yfit(d3bins)
#if only using DDD, RRR
qbins = (zetabins)/(xi1*xi2+xi2*xi3+xi3*xi1)

fig=plt.figure()
plt.plot(bins,qbins,'.')
plt.xlabel(r'$\theta$ (degrees)')
plt.ylabel('Q')
plt.title('Q vs. angle')
fig.savefig(runname+'_QvsAngle.png')
