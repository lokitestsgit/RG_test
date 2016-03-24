#!/usr/bin/python

# This code uses the logP-Mag distribution of sample A to determine M of
# the sample B, given the LogP distribution of B sample.

import astropy, os, re
from astropy.table import Table
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from astropy.table import vstack
from sklearn.neighbors.kde import KernelDensity

plt.ion()

def statuscounter(i,d):
	
	if np.remainder(i,d)==0:
		print `i`+"/"+`lenB`

def FITStablesaver(T,nfn):
	
	try:
		T.write(nfn, format='fits')
	except IOError:
		os.remove(nfn)
		T.write(nfn, format='fits')

def tableprep(fn):
	
	""" Clean the tables and add some convenient columns. """
	
	T = Table.read(fn)
	nfn = fn.replace('.fits', '_2.fits')# New table file name

	T = make_em(T)
	FITStablesaver(T, nfn)

def dataload():
	
	""" Load the calibrating sample (A) and sample with unknown 
	distances (B). """
	
	global A, B, fnA, fnB, lPcnA, lPcnB
	
	dwd = os.getcwd() # Data WD
		
	# First sample A is loaded. This is the "calibrating" sample.
	# In this case it is the OGLE III LMC small amplitude RGs.
	
	fnA = '/LMC-CalSample-cleaned_2.fits'
	A = Table.read(dwd+fnA)

	# Then sample B is loaded. For comparison/testing purposes, this is
	# again the OGLE III LMC SARGs.
	
	fnB = '/LMC-CalSample-cleaned_2.fits'
	B = Table.read(dwd+fnB)
	
	""" Fix tables so only the stars with all three good periods are 
	considered. """
	
	lPcnA = get_logPcn(A)
	lPcnB = get_logPcn(B)
	
	for cn in lPcnA:
		A = A[A[cn]>0]
	for cn in lPcnB:
		B = B[B[cn]>0]

def result_table(s, star, pdf_mag, delta_mag, deltas, KDE_mag, KDEdelta_mag, sigma, nstar):

	outTab[s]['ID'] = star['OGLEID']
	outTab[s]['WJK'] = star['ECM']
	outTab[s]['est_mag'] = pdf_mag
	outTab[s]['delta_mag'] = delta_mag
	outTab[s]['delta1'] = deltas[0]
	outTab[s]['delta2'] = deltas[1]
	outTab[s]['delta3'] = deltas[2]
	outTab[s]['KDE_mag'] = KDE_mag
	outTab[s]['KDEdelta_mag'] = KDEdelta_mag
	outTab[s]['sigma'] = sigma
	outTab[s]['nstar'] = len(C)
		
def defstuff():
	
	""" Here I define some variables and column names. """
	
	global PA, PB, col, col2, rng, xlimits, nbin, lPbw, WJK, outTab
	
	PA = ['Per1', 'Per2', 'Per3', 'Per4', 'Per5', 'Per6', 'Per7', 'Per8', 'Per9', 'Per10'] # Period columns for A sample
	PB = ['P_1', 'P_2', 'P_3'] # Period columns for B sample
	# logPB = ['logP_1', 'logP_2', 'logP_3'] 
	col = {1:'r', 2:'g', 3:'b'} 
	col2 = {1:'m', 2:'y', 3:'k'}
	rng = (8,14) # Magnitude range
	xlimits = (0.3 ,3.0) # X-axis plot limits
	bw = 0.01 # histogram bin width -- not global!
	nbin = (max(rng)-min(rng))/bw # How many bins for histogram.
	#~ lPbw = 0.025 # log period bin width
	################# CAREFUL!!!!! #####################
	lPbw = 0.025 # log period bin width
	
	outTab = Table(np.zeros((len(B), 11)), names=('ID', 'WJK', 'est_mag', 'delta_mag', 'delta1', 'delta2', 'delta3', 'KDE_mag', 'KDEdelta_mag', 'sigma', 'nstar'), dtype=('string', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64', 'float64' ))

def get_Pcn(tab):
		
	regex=re.compile("^(P_).*")
	return [m.group(0) for l in tab.colnames for m in [regex.search(l)] if m]
	
def get_logPcn(tab):
		
	regex=re.compile("^(logP_).*")
	return [m.group(0) for l in tab.colnames for m in [regex.search(l)] if m]
	
def Plimits(star, lPcn):
	
	limits = []
	
	for i in lPcn:
		lPl = star[i]*(1-lPbw)
		lPh = star[i]*(1+lPbw)
		limits.append((lPl, lPh))

	return limits

def make_em(T): 
	# make a replacement for the extinction-corrected magnitude column. T short for temptab
	
	ecm = Table.Column(T['Kf'], name='ECM', description='Kf magnitude -- the best K band mag')
	
	#~ ecm = Table.Column(T['Kf']-0.686*(T['Jf']-T['Kf']), name='ECM', description='WJK') # 2MASS extinction-corrected magnitude column, WJK in this case
			
	T.add_column(ecm)
	return T
	 
	ecm = Table.Column(ecm, name='ECM')
	T.add_column(ecm)
	return T

def PLdeco1():
	
	plt.ylim(max(rng), min(rng))
	plt.xlim(xlimits)
	plt.xlabel("$log P$", fontsize=16)
	plt.ylabel("$K_{S}$", fontsize=16)
	plt.tick_params(labelsize=16)

	plt.draw()
	plt.show()

def PLdeco3():
	
	plt.xlabel("$K_{S}$", fontsize=16)
	plt.ylabel("$p$", fontsize=16)
	plt.draw()
	plt.show()

def PLplot(limits, star, bin_edges, pdf, pdf_mag, pdf_max, X_plot, log_dens_exp, KDE_mag, sigma):
	
	lPcnA = get_logPcn(A)
	lPcnB = get_logPcn(B)
	lPcnC = get_logPcn(C)
	
	fig = plt.figure("PL diagram")
	plt.clf()
	ax = fig.add_subplot(1,3,1)
	plt.tick_params(labelsize=16)

	for j in lPcnA:
		plt.scatter(A[j], A['ECM'], c='k', s=0.05, zorder=1, alpha=0.2)

	for i, j in enumerate(limits):
		plt.vlines(j[0], rng[0], rng[1], colors=col[i+1])
		plt.vlines(j[1], rng[0], rng[1], colors=col[i+1])
		ax.add_patch( patches.Rectangle( (j[0], min(rng)), j[1]-j[0], np.ptp(rng), color=col[i+1], alpha=0.25 ) )

	for j in lPcnB:
		plt.plot(star[j], star['ECM'], marker='x', mec='yellow', mew=2, zorder=4, alpha=1)
			
	plt.ylim(max(rng), min(rng))
	plt.xlim(xlimits)
	plt.xlabel("$log P$", fontsize=16)
	plt.ylabel("$K_{S}$", fontsize=16)

	ax = fig.add_subplot(1,3,2)

	for j in lPcnC:
		plt.scatter(C[j], C['ECM'], c='k', s=1, zorder=1, alpha=0.7)
		
	for i, j in enumerate(limits):
		plt.vlines(j[0], rng[0], rng[1], colors=col[i+1])
		plt.vlines(j[1], rng[0], rng[1], colors=col[i+1])
		ax.add_patch( patches.Rectangle( (j[0], min(rng)), j[1]-j[0], np.ptp(rng), color=col[i+1], alpha=0.25 ) )
	plt.hlines(star['ECM'], xlimits[0], xlimits[1], linestyles='dashed', color='y', lw=2)
	plt.hlines(pdf_mag, xlimits[0], xlimits[1], linestyles='dashed', color='c', lw=2)
	
	plt.ylim(max(rng), min(rng))
	plt.xlim(xlimits)
	plt.xlabel("$log P$", fontsize=16)
	plt.ylabel("$K_{S}$", fontsize=16)
	plt.tick_params(labelsize=16)

	
	ax = fig.add_subplot(1,3,3)
	
	plt.plot(bin_edges, np.hstack((pdf[0],pdf)), drawstyle='steps', zorder=2, color='k', alpha=1, lw=2, linestyle='-')
	plt.vlines(star['ECM'], 0, pdf_max, linestyles='-', color='y', lw=2)
	plt.vlines(pdf_mag, 0, pdf_max, linestyles='-', color='c', lw=2)
	ax.add_patch( patches.Rectangle( (min((star['ECM'], pdf_mag)), 0), abs(star['ECM']-pdf_mag), pdf_max, color=col[i+1], alpha=0.2, zorder=1 ) )
	
	plt.xlabel("$K_{S}$", fontsize=16)
	plt.ylabel("$p$", fontsize=16)
	plt.tick_params(labelsize=16)

	ax2=plt.twinx()

	for tick in ax2.yaxis.get_major_ticks():
		tick.label1On = False
		tick.label2On = True
	
	plt.plot(X_plot[:,0], log_dens_exp, '-', color='r', lw=2, label='KDE', alpha=0.75)
	plt.ylabel('Normalized density',fontsize=16)
	plt.tick_params(labelsize=16)
	ax2.yaxis.set_label_position("right")
	
	plt.draw()
	plt.show()

def PLplotT1(i, limits, star, bin_edges, pdf, pdf_mag, pdf_max):

	global fig	
	lPcnC = get_logPcn(C)
	fig = plt.figure("PL diagram")
	plt.clf()
	ax = fig.add_subplot(1,3,1)

	plt.scatter(jA['lP'], jA['ECM'], c='k', s=0.05, zorder=1, alpha=0.2)

	for j, l in enumerate(limits):
		plt.vlines(l[0], rng[0], rng[1], colors=col[j+1])
		plt.vlines(l[1], rng[0], rng[1], colors=col[j+1])
		ax.add_patch( patches.Rectangle( (l[0], min(rng)), l[1]-l[0], np.ptp(rng), color=col[j+1], alpha=0.25 ) )

	for j in lPcnB:
		plt.plot(star[j], star['ECM'], marker='x', mec='yellow', mew=2, zorder=4, alpha=1)

	PLdeco1()
	
	ax = fig.add_subplot(1,3,2)

	#~ for j in lPcnC:
		#~ plt.scatter(C[j], C['ECM'], c='k', s=1, zorder=1, alpha=0.7)
	plt.scatter(C['lP'], C['ECM'], c='k', s=0.1, zorder=1, alpha=0.25)
		
	for j, l in enumerate(limits):
		plt.vlines(l[0], rng[0], rng[1], colors=col[j+1])
		plt.vlines(l[1], rng[0], rng[1], colors=col[j+1])
		ax.add_patch( patches.Rectangle( (l[0], min(rng)), l[1]-l[0], np.ptp(rng), color=col[j+1], alpha=0.25 ) )
	plt.hlines(star['ECM'], xlimits[0], xlimits[1], linestyles='dashed', color='y', lw=2)
	plt.hlines(pdf_mag, xlimits[0], xlimits[1], linestyles='dashed', color=col[i+1], lw=2)
	
	PLdeco1()
	
	ax = fig.add_subplot(1,3,3)
	
	plt.plot(bin_edges, np.hstack((pdf[0],pdf)), drawstyle='steps', zorder=2, color=col[i+1], alpha=1, lw=2, linestyle='-')
	plt.vlines(pdf_mag, 0, pdf_max, linestyles='-', color=col[i+1], lw=2)
	
	PLdeco3

	plt.draw()
	plt.show()
	
def PLplotT2(i, limits, star, bin_edges, pdf, pdf_mag, pdf_max):

	ax = fig.add_subplot(1,3,2)

	plt.scatter(C['lP'], C['ECM'], c='k', s=0.1, zorder=1, alpha=0.25)
		
	plt.hlines(pdf_mag, xlimits[0], xlimits[1], linestyles='dashed', color=col[i+1], lw=2)
	
	PLdeco1()
	
	ax = fig.add_subplot(1,3,3)
	
	plt.plot(bin_edges, np.hstack((pdf[0],pdf)), drawstyle='steps', zorder=2, color=col[i+1], alpha=1, lw=2, linestyle='-')
	plt.vlines(pdf_mag, 0, pdf_max, linestyles='-', color=col[i+1], lw=2)
	
	PLdeco3()

	plt.draw()
	plt.show()

def PLplotT3(i, limits, star, bin_edges, pdf, pdf_mag, pdf_max):
	
	ax = fig.add_subplot(1,3,2)
	plt.scatter(C['lP'], C['ECM'], c='k', s=0.1, zorder=1, alpha=0.25)
	plt.hlines(pdf_mag, xlimits[0], xlimits[1], linestyles='dashed', color='k', lw=2)
	
	PLdeco1()
	
	ax = fig.add_subplot(1,3,3)
	
	plt.plot(bin_edges, np.hstack((pdf[0],pdf)), drawstyle='steps', zorder=2, color='k', alpha=1, lw=2, linestyle='-')
	plt.vlines(star['ECM'], 0, pdf_max, linestyles='-', color='y', lw=2)
	plt.vlines(pdf_mag, 0, 10*pdf_max, linestyles='-', color='k', lw=2)
	ax.add_patch( patches.Rectangle( (min((star['ECM'], pdf_mag)), 0), abs(star['ECM']-pdf_mag), 10*pdf_max, color=col[i+1], alpha=0.2, zorder=1 ) )
	
	PLdeco3()

	plt.draw()
	plt.show()

def find_sim_stars(limits, star):

	""" Find stars from the calibrating sample that have all three periods
	matching the periods of the test star (within 0.05*log10(P)).  """

	global C
	
	st = np.zeros((len(A),len(lPcnB)))
	
	for l in limits:
		for k, cn in enumerate(lPcnB):
			st[:,k][(A[cn]>l[0]) & (A[cn]<l[1])]=1
	
	C = A[(st[:,0]>0) & (st[:,1]>0) & (st[:,2]>0)] # C is the table with the similar stars -- for Laurent's version
	C = C[C['OGLEID']!=star['OGLEID']] # Remove the "star" in case of self-reconstruction

	return C
		
def makehist(col):
	# Make a pdf from the selected stars. Count the number of stars in
	# 0.01 mag wide bins and normalize.
	
	hist_val, bin_edges = np.histogram(col, bins=nbin, range=rng, density=True) # density=True: area under the curve normalized to unity.
	hvsn = np.convolve(hist_val, np.ones((10))/10, mode='same') # Smooth the histogram
	return bin_edges, hvsn

def makeKDE(m):

	m = m[(m>-100) & (m<100)]
	m = m[:, np.newaxis] # Training data
	l = len(m)
	sigma = np.std(m)
	kdebw = (1.*4/3*sigma**5/ l)**(1./5.)
	
	try:
		X_plot = np.linspace(rng[0], rng[1], 1000)[:, np.newaxis]
		kde = KernelDensity(kernel='gaussian', bandwidth=kdebw).fit(m)
		log_dens = kde.score_samples(X_plot)
		log_dens_exp = np.exp(log_dens)
		KDE_mag = np.float(X_plot[np.argmax(log_dens_exp)])
	except ValueError:
		log_dens_exp = np.ones(len(X_plot[:,0]))*-99.99
		KDE_mag, sigma = -99.99, -99.99
	return X_plot, log_dens_exp, KDE_mag, sigma
	
def mkjoinTab(A):
	
	global jA
	for i, cn in enumerate(lPcnA):
		if i==0:
			jA = A['ECM', cn]
			jA.rename_column(cn, 'lP')
		else:
			t2 = A['ECM', cn]
			t2.rename_column(cn, 'lP')
			jA = vstack([jA, t2])

	return jA

def diagnostics(pdf, bin_edges, star):
	
	pdf_mag = bin_edges[np.argmax(pdf)]
	pdf_max = max(pdf)
	tmag = star['ECM'] # "True" mag
	delta_mag = pdf_mag-tmag

	return pdf_mag, pdf_max, delta_mag

def KDEdiagnostics(KDE_mag, star, sigma):
	
	tmag = star['ECM'] # "True" mag
	KDEdelta_mag = KDE_mag-tmag
	print "True mag:", tmag, "KDE Test mag:", KDE_mag, "delta mag:", KDEdelta_mag, "sigma:", round(sigma,2)
	return KDEdelta_mag

def toploop():
	""" Stars with *all* the periods matching those of the test star 
	will be selected. """

	global lenB, deltas
	lenB = len(B)
	jA = mkjoinTab(A)
	
	for s, star in enumerate(B[5000:]):
		deltas = [-99.99,-99.99,-99.99]
		limits = Plimits(star, lPcnB)

		# Restrictive version
		C = find_sim_stars(limits, star)
		if len(C)>=10:
			bin_edges, pdf = makehist(C['ECM'])
			pdf_mag, pdf_max, delta_mag = diagnostics(pdf, bin_edges, star)
			X_plot, log_dens_exp, KDE_mag, sigma = makeKDE(C['ECM'])
			KDEdelta_mag = KDEdiagnostics(KDE_mag, star, sigma)
			PLplot(limits, star, bin_edges, pdf, pdf_mag, pdf_max, X_plot, log_dens_exp, KDE_mag, sigma)
			result_table(s, star, pdf_mag, delta_mag, deltas, KDE_mag, KDEdelta_mag, sigma, len(C))
			raw_input()
		else:
			result_table(s, star, -999.99, -999.99, (-999.99, -999.99, -999.99), -999.99, -999.99, -999.99, len(C))
		statuscounter(s, 100)
		#~ raw_input()

dataload()
defstuff()
toploop()
#FITStablesaver(outTab, "BLG_estM-wEs.fits")
