# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 20:04:17 2015

pp_tools.py
Code for the Point Process Toolbox (https://github.com/gmf05/pp_tools)
Ported from MATLAB

@author: gmf
"""

""" MODULES """
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.stats
from statsmodels.distributions.empirical_distribution import ECDF
#plt.style.use('ggplot')
#plt.style.use('fivethirtyeight')
plt.style.use('bmh')

""" GLOBAL VARAIBLES """
PLOT_COLOR = 'b'
DO_CONF_INT = False
DO_MASK = False
RESCALING = 'exp'
WINDOW = np.array([2, 0.2])
PRECISION = 3 # Number of decimal points to display

""" POINT PROCESS DATA """
class data:
  name = ''
  labels = []
  dn = []  
  time = []
  len_time = []
  dt = []
  sampling_rate = []
  num_channels = 0
  
  def __init__(self, dn, time):
    self.dn = np.array(np.matrix(dn))
    self.num_channels = self.dn.shape[0]
    self.time = np.array(time)
    self.len_time = len(time)
    self.dt = (time[1] - time[0]) * 1.0 # cast to float
    self.sampling_rate = 1.0 / self.dt
    #Print 'Made point process data object'
  
  def __repr__(self):
    return self.__str__()
    
  def __str__(self):
    print ''
    print 'Point Process Data'
    print 'Name : %s' % self.name    
    print 'Number channels : %d' % self.num_channels
    print 'Start time : %s ' % str(self.time[0])
    print 'End time : %s' % str(self.time[-1])
    print 'Number time points : %d' % self.len_time
    print 'Sampling Rate: %f Hz' % self.sampling_rate
    print 'Time step dt : %s' % self.dt
    return ''
    
  def append(self, dn, axis=0):
    self.dn = np.vstack((self.dn, dn))
    self.num_channels = len(self.dn)
     
  def sub_time(self, t_start, t_end):
    junk,i_on = find_nearest(self.time, t_start)
    junk,i_off = find_nearest(self.time, t_end)
    ti = range(i_on, i_off)
    return data(self.dn[:,ti], self.time[ti])
  
  #def sub_time(self, ti): 
  #  return data(self.dn[:,ti], self.time[ti])  
  
  def sub_data(self, channel_list):
    return data(self.dn[channel_list,:],self.time)

  def downsample(self, downsample_factor):
    print 'Downsampling point process data by a factor of %d...' % downsample_factor
    self.sampling_rate /= downsample_factor
    self.dt *= downsample_factor
    self.time = np.array([self.time[i] for i in range(0,self.len_time,downsample_factor)])
    self.len_time = len(self.time)
    self.dn = cumdownsample(self.dn, downsample_factor)
    if np.any(self.dn>1): print 'Warning! Downsampling has caused values > 1.'
    return self
    
  def plot(self, plot_type='raster'):

    #plt.figure() # Open new figure
    
    # Partition time axis into windows
    win_bins = np.round( WINDOW[0] / self.dt )
    time_window = np.arange(self.time[0], self.time[-1]-WINDOW[0], WINDOW[1])
    num_windows = len(time_window)
    
    if plot_type=='raster':
      for i in range(self.num_channels):
        spikei = np.nonzero(self.dn[i,:])[0] # Is this wrong??
        #spikei = np.nonzero(self.dn[i,:])[1] # Want non-zero COLUMNS
        #plt.plot(self.time[spikei], i*np.ones(len(spikei)), PLOT_COLOR + '.')
        plt.plot(self.time[spikei], i*np.ones(len(spikei)), PLOT_COLOR + '|', markersize=24)
      plt.title('%s Raster Plot' % self.name)
      plt.xlabel('Time [s]')
      plt.ylabel('Channel')
      plt.gca().invert_yaxis()

    elif plot_type=='raster-marks':
      0 # TODO: ADD
      
    elif plot_type=='rate':
      rate_window = np.zeros([self.num_channels, num_windows])
      count=0
      for n in range(num_windows):
        rate_window[:, n] = np.sum(self.dn[:, count:count+win_bins], axis=1)
        count += win_bins
      print time_window.shape
      print rate_window.shape
      #plt.plot(np.matrix(time_window), rate_window)
      plt.plot(time_window, rate_window.T)
      plt.title('%s Firing Rate' % self.name)
      plt.xlabel('Time [s]')
      plt.ylabel('Spikes / sec')
      
    elif plot_type=='psth':
      plt.plot(self.time, np.sum(self.dn, axis=0))
      plt.title('%s PSTH' % self.name)
      plt.xlabel('Time [s]')
      plt.ylabel('Number of spikes')
      
    elif plot_type=='heat':
      2 # TODO: ADD
      
    elif plot_type=='isi':
      3 # TODO: ADD
      
    elif plot_type=='isi-heat':
      3 # TODO: ADD
      
    else:
      raise ValueError('Unknown plot type %s' % plot_type)
    
""" POINT PROCESS MODEL """
class stats:
  dfe = 0
  s = 1
  estdisp = 0
  se = []
  t = []
  p = []
  # Add residuals??
  weights = []
    
  def __str__(self):
    print 'Point process model statistics'
    print 'Degrees of freedom for error : %s' % str(self.dfe)
    print 'Standard Error : %s' % str(np.around(self.se, decimals=PRECISION))
    print 't-statistics for coefficients : %s' % str(np.around(self.t, decimals=PRECISION))
    print 'p-values for coefficients : %s' % str(np.around(self.p, decimals=PRECISION))
    # Add residuals??
    return ''
    
  def __repr__(self):
    return self.__str__()
    
  def isempty(self): return bool(self.dfe==0)


class model:
  name = ''
  coeff = []
  covar = []
  design = []
  response = []
  intensity = []
  time = []
  fit_method = 'irls'
  link = 'log'
  params = []
  stats = stats()
  rescaled_isi = []
  log_likelihood = np.nan
  deviance = np.nan
  aic = np.nan
  ks = []
  
  def __init__(self):
    print 'Made pp_model object'
    
  def __repr__(self):
    return self.__str__()
    
  def __str__(self):
    print ''
    print 'Point Process Model'
    print 'Name : %s' % self.name
    print 'Coefficients : %s' % str(np.around(np.squeeze(np.array(self.coeff)), decimals=PRECISION))
    print 'Link function : %s' % self.link
    print 'Fit method : %s' % self.fit_method
    print 'Number spikes %d' % np.sum(self.response)
    if ~self.stats.isempty(): print 'p-values : %s' % str(np.around(self.stats.p, decimals=PRECISION))
    print 'Log Likelihood : %s' % str(self.log_likelihood)
    print 'Deviance : %s' % str(self.deviance)
    print 'Akaike Information Criterion AIC : %s' % str(self.aic)
    try:
      print 'KS statistic : %s' % str(self.ks[0])
      print 'KS statistic bound (95% conf int) : %s ' % str(self.ks[1])
    except:
      0
    return ''

  def fit(self, d, p): # d = pp.data, p = pp.params

    print 'Fitting pp model...'
    
    # Make design matrix
    print 'Making design matrix...'
    self.design, self.response = self.makeXY(d, p)
    
    # Estimate model parameters according to fit_method
    if self.fit_method=='irls':
      print 'Running IRLS...'
      self.coeff,self.covar,self.stats = irls(self.design, self.response, self.link)
    else:
      raise ValueError('Unknown fit method')
    
    # NOTE: NEED TO CHANGE 'EXP' to inverse link in general...
    self.intensity = np.exp( np.dot(self.design , self.coeff ) )
    self.time = d.time
    self.params = p
    #self.calc_gof()
    return self
    
  def makeXY(self, d, p):
    
    # Initialize stuff
    y = d.dn[p.response,:].T  # data
    T = d.len_time # number of time steps
    N = p.num_covar() # number of covariates
    X = np.zeros([T,N]) # design matrix
    burnin = p.get_burnin()
    
    for i in range(p.num_covar_groups()):

      ind = p.index[i]
      knots = p.covariates[i].knots
      nknots = len(knots)
      count = 0 # used to keep track of index in time (i.e. which row)
      
      # INTERCEPT / BASELINE SPIKE PROBABILITY TERM(S)
      if p.covariates[i].channel<0:
        bins_per_knot = np.round(np.diff(knots) * T)
        if p.covariates[i].basis=='indicator':
          for n in range(nknots-1):
            X[count:count+bins_per_knot[n], ind[n]] = 1
            count += bins_per_knot[n]
        elif p.covariates[i].basis=='spline':
          s_coeff = p.tension_matrix()
          for n in range(nknots-1):
            nt = bins_per_knot[n]
            alphas = 1/(nt-1) * np.arange(0,nt,1)
            X[count:count+nt, ind[n]:ind[n]+4] = np.matrix([alphas**3, alphas**2, alphas, np.ones([nt,1])]).T * s_coeff
            count += nt
        else:
          raise ValueError('Unknown basis function %s' % p.covariates[i].basis)
      # PROBABILITIES RELATED TO SPIKE HISTORIES
      # e.g. self-history (rhythms)
      # or ensemble-history (network coupling)
      else:
        # Compute sum of data if number of channels > 1
        #if len(p.covariates[i].channel)>1:
        if False:
          di = np.sum(d.dn[p.covariates[i].channel,:], axis=0)
        else:
          di = d.dn[p.covariates[i].channel,:]
       
        if p.covariates[i].basis=='indicator':
          for n in range(nknots):
            X[burnin:, ind[n]] =  di[burnin-knots[n]:]
            
        elif p.covariates[i].basis=='spline':
          Xs = p.splineX(i)
          onset = p.covariates[i].knots[0]
          offset = p.covariates[i].knots[-1]+1
          spikei = d.dn[p.covariates[i].channel,:].nonzero()[0]
          for s in spikei:
            bins_to_end = T - s
            if s + offset > T:
              X[s+onset : s+bins_to_end, ind] += di[s] * Xs[:bins_to_end-onset, :]              
            else:              
              X[s+onset : s+offset, ind] += di[s] * Xs
        else:
          raise ValueError('Unknown basis function %s' % p.covariates[i].basis)

    # Trim burn-in period, save X & y
    X = X[burnin:, :]
    y = y[burnin:]

    #self.design = X
    #self.response = y
    return X,y

  def plot(self):

    #plt.figure()
    
    num_sec = self.time[-1] - self.time[0]
    dt = self.time[1] - self.time[0]
    bins_per_ms = 1/(dt*1e3)
    Nc = len(self.params.covariates)
    for c in range(Nc):
      knots = self.params.covariates[c].knots
      plt.subplot(Nc,1,c+1)
      plt.title(self.params.covariates[c].name)
      ci = self.params.index[c]
      if self.params.covariates[c].basis=='indicator':
        if self.params.covariates[c].channel<0:
         for n in range(len(knots)-1):
           #plt.plot(knots[[n,n+1]], np.exp( self.coeff[p.index[c][n],0] ) / dt * np.ones(2), PLOT_COLOR)
           #print self.coeff[p.index[c][n],0]
           plt.plot(knots[[n,n+1]] * num_sec + self.time[0], np.exp( self.coeff[self.params.index[c][n],0] ) / dt * np.ones(2), PLOT_COLOR)
         plt.xlabel('Time [s]')
         plt.ylabel('[Hz]')   
        else:
          plt.plot(knots / bins_per_ms, self.coeff[ci])
          # compute,plot upper yhi
          # compute,plot lower ylo
          plt.xlabel('Lag Time [ms]')
          plt.ylabel('Mod.')   
        if DO_CONF_INT:
          0          
          #yhi = X0*coeff + Z*np.matrix(np.sqrt(np.diag(np.dot(X0,np.dot(cov,X0.T))))).T # Z=2: 95% CI
          #ylo = X0*coeff - Z*np.matrix(np.sqrt(np.diag(np.dot(X0,np.dot(cov,X0.T))))).T # Z=2 : 95% CI
          #plot
      elif self.params.covariates[c].basis=='spline':

        if DO_CONF_INT:
            tspline,y,yhi,ylo = plotspline(knots, self.coeff[ci], self.covar[np.ix_(ci,ci)])
        else:
            tspline,y = plotspline(knots, self.coeff[ci])
        
        if self.params.covariates[c].channel<0:          
          tspline = tspline * num_sec + self.time[0]
          scale_factor = 1.0 / dt
          plt.xlabel('Time [s]')
          plt.ylabel('[Hz]')   
        else:
          tspline = tspline / bins_per_ms
          scale_factor = 1
          plt.xlabel('Lag Time [ms]')
          plt.ylabel('Mod.')
        y = np.exp(y) * scale_factor
        plt.plot(tspline, y, PLOT_COLOR)
        if DO_CONF_INT:
          yhi = np.exp(yhi) * scale_factor
          ylo = np.exp(ylo) * scale_factor
          #############
          plt.plot(tspline,yhi,'r--')

          plt.plot(tspline,ylo,'r--')
          #plt.fill_between(tspline, ylo, yhi, alpha=0.5, edgecolor=PLOT_COLOR, facecolor=PLOT_COLOR)          
          #############
      else:
        raise ValueError('Basis function %s not recognized' % self.params.covariates[c].basis)
  
  ## return intensity + bounds on intensity?
  def val(self):    
    return ''
  ## return intensity + bounds on intensity?
    
  def calc_gof(self):
    if self.link == 'log':
      def link_pdf(y1, y2): return scipy.stats.poisson.pmf(y1, y2)
    elif self.link == 'logit':
      def link_pdf(y1, y2): return scipy.stats.binom.pmf(y1, y2)
    else:
      raise ValueError('Unrecognized link function %s' % self.link)

    #self.log_likelihood=0
    #self.deviance=0
    
    #self.log_likelihood = np.sum(np.log(link_pdf(self.response, np.squeeze(np.array(self.intensity)))))
    self.log_likelihood = 0
    self.deviance = 2*(np.sum(np.log(link_pdf(self.response,self.response))) - self.log_likelihood)

    self.aic = self.deviance + 2 * len(self.coeff)
    self.rescaled_isi = rescale_isi(self.response, self.intensity)
    self.ks = kstest(self.response, self.intensity)
    return self        
    
  def gof_plot(self):

    plt.figure() # Open new figure
    
    # Conditional Intensity Function
    T = len(self.intensity)
    plt.subplot(131)
    plt.plot(self.time[-T:], self.intensity, PLOT_COLOR)
    plt.title('Conditional Intensity Function')
    plt.xlabel('Time [s]')
    plt.ylabel('Spikes / bin')
       
    # Histogram of Rescaled ISIs
    plt.subplot(132)
    y,x=np.histogram(self.rescaled_isi,bins=np.arange(0,1.05,0.1))
    dx = x[1]-x[0]
    plt.bar(x[:-1],y, width=dx,color=PLOT_COLOR)
    plt.title('Histogram of Rescaled ISIs')
    plt.xlabel('Inter-Spike Interval')
    plt.ylabel('Count')
    
    # Kolmogorov-Smirnov (KS) Plot
    plt.subplot(133)
    self.ks_plot()
    
  def residual_plot(self):
    cum_cif = np.cumsum(np.array(self.intensity))
    cum_spikes = np.cumsum(np.array(self.response))
    resd = cum_spikes - cum_cif
    plt.plot(self.time[:len(resd)], resd, PLOT_COLOR)
    plt.title('Residual Plot')
    plt.xlabel('Time [s]')
    plt.ylabel('Observed Spikes - Expected Spikes')
    
  def qq_plot(self):
    z = np.sort(self.rescaled_isi)
    if RESCALING=='exp':
      def rsinv(x): return scipy.stats.uniform.ppf(x)
    else:
      def rsinv(x): return scipy.stats.expon.ppf(x, 1)
    I = np.arange(0,1,0.05)
    qx = rsinv(I)
    qz = np.array([np.percentile(z, i) for i in 20*I])
    print len(qx)
    print len(qz)
    plt.plot(qx, qz, PLOT_COLOR)
    plt.plot(I,I,'r--')
    #ci=0 # ??    
    #plt.plot(I,I+ci,'r--')
    #plt.plot(I,I-ci,'r--')
    plt.title('Quantile-Quantile (QQ) Plot')
    plt.xlabel('Empirical Quantiles')
    plt.ylabel('Theoretical Quantiles')

  def autocorr_plot(self):
    ac,lags = autocorr(self.rescaled_isi)
    plt.plot(lags, ac, PLOT_COLOR)
    plt.xlabel('Number of spikes prev')
    plt.ylabel('Autocorrelation')

  def ks_plot(self):
    z = rescale_isi(self.response, self.intensity)
    if RESCALING=='exp':
      def rs_cdf(x): return scipy.stats.uniform.cdf(x)
    else:
      def rs_cdf(x): return scipy.stats.expon.cdf(x, 1)
    ecdf = ECDF(np.sort(z))
    acdf = rs_cdf(ecdf.x)
    #ks_stat = np.round( np.max(np.abs(acdf - ecdf.y)) , 3)
    plt.plot(ecdf.y, acdf, PLOT_COLOR) # empirical vs actual cdf
    ks_ci = np.round( 1.36/np.sqrt(len(z)), 3 )
    xx = np.arange(0,1,0.05)
    plt.plot(xx, xx, 'r--') # diagonal : "perfect" model
    plt.plot(xx, xx + ks_ci, 'r--') # confidence bounds
    plt.plot(xx, xx - ks_ci, 'r--') # confidence bounds
    plt.title('Kolmogorov-Smirnov (KS) Plot')
    plt.xlabel('Empirical CDF')
    plt.ylabel('Theoretical CDF')
    #return ecdf,acdf

""" POINT PROCESS PARAMETERS """
class params:
  covariates = []
  index = []
  response = 0
  fit_method = 'irls'
  noise = []
  window = []
  downsample_est = 1
  
  def __repr__(self):
    return self.__str__()
    
  def __str__(self):
    #print 'Coefficients'
    #
    print ''
    print 'Point Process Parameters'
    print 'Covariate Names : %s' % ', '.join([self.covariates[i].name for i in range(self.num_covar_groups())])
    print 'Covariate Knots : %s' % ', '.join([str(self.covariates[i].knots) for i in range(self.num_covar_groups())])
    print 'Covariate Channels : %s' % ', '.join([str(self.covariates[i].channel) for i in range(self.num_covar_groups())])
    print 'Covariate Basis : %s' % ', '.join([self.covariates[i].basis for i in range(self.num_covar_groups())])
    print 'Fit method : %s' % self.fit_method
    print 'Number covariates : %s ' % str(self.num_covar())
    print 'Response channel : %s ' % str(self.response)

    try:
      print 'KS statistic : %s' % str(self.ks[0])
      print 'KS statistic bound (95% conf int) : %s ' % str(self.ks[1])
    except:
      0
    return ''
  
  def add_covar(self, name, channel, knots, basis):
    if len(self.index)==0:
      istart = 0
    else:
      istart = self.index[-1][-1] + 1
    istop = istart + len(knots) + 2 * (basis=='spline') - 1*(channel<0)*(basis=='indicator')
    self.index.append( range(istart,istop) )
    self.covariates.append( covariate( name, channel, knots, basis ) )
    
  def num_covar(self):
    return self.index[-1][-1]+1 # Seems we *DO* want +1 here
  
  def num_covar_groups(self):
    return len(self.covariates)
    
  def get_burnin(self):
    burnin = 0
    for c in self.covariates:
      if c.channel>=0: burnin = np.max([burnin, c.knots[-1]])
    return int(burnin)
    
  def tension_matrix(self, s=0.5):
    return np.array([[-s, 2.0-s, s-2.0, s],[2.0*s, s-3.0, 3.0-2.0*s, -s],[-s, 0, s, 0],[0, 1, 0, 0]])
  
  def splineX(self, i):
    knots = np.copy( self.covariates[i].knots )
    nknots = len(knots)  
    nlags = knots[-1] - knots[0] + 1
    if  nlags>1000: dtau = 1
    elif nlags>1: dtau = 1
    else: dtau = 1e-2
    s_coeff = self.tension_matrix()
    tau = np.arange(knots[0], knots[-1]+dtau, dtau)
    ntau = len(tau)
    X0 = np.zeros([ntau,nknots+2])    
    count=0
    for i in range(nknots-1):
      alphas = np.linspace(0,1,num=knots[i+1]-knots[i])
      Na = len(alphas)
      A = np.matrix([alphas**3, alphas**2, alphas, np.ones([Na,1])]).T * s_coeff
      X0[count:count+Na, i:i+4] = A
      count += Na
    X0[-1, i:i+4] = np.matrix([1,1,1,1]) * s_coeff
    return X0
  
  def __init__(self):
    self.covariates = []
    self.index = []
    self.response = 0
    self.fit_method = 'irls'
    self.noise = []
    self.window = []
    self.downsample_est = 1
    print 'Made pp_params object'

""" MODEL COVARIATE """
class covariate:
  name = ''
  channel = 0
  knots = []
  basis = ''
  def __init__(self, name, channel, knots, basis):
    self.name = name
    self.channel = channel
    self.knots = np.array(knots)
    self.basis = basis

""" HELPER FUNCITONS """   
def autocorr(x):
  x = x - np.mean(x) # make unbiased  
  nlags = len(x)
  lags = np.arange(0,nlags)
  ac = np.array( [1]+[np.corrcoef(x[:-i], x[i:])[0,1] for i in range(1, nlags)] )
  return ac,lags
  
def cumdownsample(x, downsample_factor):  
  Nr,Nc = x.shape
  NT = len(range(0,Nc,downsample_factor))
  y = np.zeros([Nr,NT])
  for n in range(NT-1):
    y[:,n] = np.sum(x[:, (n-1)*downsample_factor:n*downsample_factor-1], axis=1)
  return y

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx],idx

def plotspline(knots, coeff, cov=0, s=0.5, Z=2):
  nknots = len(knots)  
  #nlags = knots[-1] - knots[0] + 1
  #if  nlags>1000: dtau = 1
  #elif nlags>1: dtau = 1
  #else: dtau = 1e-2
  if knots[-1]-knots[0]<5: dtau = 1e-3
  else: dtau = 1
  s_coeff = np.matrix([[-s, 2.0-s, s-2.0, s],[2.0*s, s-3.0, 3.0-2.0*s, -s],[-s, 0, s, 0],[0, 1, 0, 0]])
  tau = np.arange(knots[0], knots[-1]+dtau, dtau)
  ntau = len(tau)
  X0 = np.zeros([ntau,nknots+2])
  
  count=0
  for i in range(nknots-1):
    Na = np.round((knots[i+1]-knots[i])/dtau)
    alphas = np.linspace(0,1,num=Na)
    A = np.dot(np.matrix([alphas**3, alphas**2, alphas, np.ones([Na,1])]).T , s_coeff)
    X0[count:count+Na, i:i+4] = A
    count += Na
  X0[-1, i:i+4] = np.dot(np.matrix([1,1,1,1]) , s_coeff)
  y = np.dot(X0, coeff)
  ## COMPUTE CONF INTERVALS! (IF DESIRED)
  if DO_CONF_INT:
    yhi = np.dot(X0,coeff) + Z*np.matrix(np.sqrt(np.diag(np.dot(X0,np.dot(cov,X0.T))))).T # Z=2: 95% CI
    ylo = np.dot(X0,coeff) - Z*np.matrix(np.sqrt(np.diag(np.dot(X0,np.dot(cov,X0.T))))).T # Z=2 : 95% CI
    if DO_MASK: 
      i1,j1 = np.where(ylo<0)
      i2,j2 = np.where(yhi>0)
      i3 = np.intersect1d(np.array(i1),np.array(i2))
      y[i3] = 0
  
  if DO_CONF_INT: return tau,y,yhi,ylo
  else: return tau,y

def irls(X,y,link):  
  #
  # CAN WE FORCE COLUMN VECTOR IN BETTER WAY??? VVV
  y = np.matrix(y).T # Want column vector!!
  # CAN WE FORCE COLUMN VECTOR IN BETTER WAY??? ^^^
  #
  if link=='identity':
    def link_fun(x): return x
    def ilink_fun(x): return x
    def ilink_prime(x): return np.ones_like(x)
    def sqrtvar_fun(x): return np.ones_like(x)
    mu = np.copy(y)  
  elif link=='log':
    def link_fun(x): return np.log(x)
    def ilink_fun(x): return np.exp(x)
    def ilink_prime(x): return np.divide(1.0, x)
    def sqrtvar_fun(x): return np.sqrt(x)
    mu = np.copy(y) + 0.25
  elif link=='logit': # TODO: NEEDS UPDATE!!!
    def link_fun(x): return x
    def ilink_fun(x): return x
    def ilink_prime(x): return np.ones_like(x)
    def sqrtvar_fun(x): return np.ones_like(x)
    mu = np.copy(y)  
  else:
    raise ValueError('Link function `%s` not recognized' % link)

  # Initialize arrays    
  NT,Np = X.shape
  pwts = np.ones([NT,1])
  b = np.ones([Np,1])
  C = np.eye(Np)
  R = np.eye(Np)
  eta = ilink_fun(mu)
  
  # Stopping parameters
  epsilon = 1e-6
  offset = 1e-3
  iter_lim = 100
  
  for iter in range(iter_lim):    
    #print 'IRLS Iteration %d' % iter
    #z = eta - offset + (y-mu) * ilink_prime(mu)   
    z = eta - offset + np.multiply((y-mu), ilink_prime(mu))
    b_old = b.copy()
    R_old = R.copy()
    deta = ilink_prime(mu)
    sqrtirls = np.multiply( np.abs(deta), sqrtvar_fun(mu) )
    sqrtw = np.divide(np.sqrt(pwts), sqrtirls)
    
    # Orthogonal (QR) decomposition of Xw
    # Avoids forming the product Xw.T * Xw
    #zw = z * sqrtw
    #Xw = X * sqrtw[:, [0 for i in range(Np)]]
    zw = np.multiply(z, sqrtw)    
    Xw = np.multiply(X, sqrtw[:, [0 for i in range(Np)]])
    Q,R = np.linalg.qr(Xw)
    try:
      b,Z1,Z2,Z3 = np.linalg.lstsq(R, np.dot(Q.T, zw))
    except ValueError:
      b = b_old.copy()
      R = R_old.copy()
      print 'ValueError at iteration %d' % iter
      break
    
    # Checks:
    # 1. Have we encountered numerical error/flat likelihood?
    cond_num = np.linalg.cond(R)
    if cond_num < 1e-8 or np.any(np.isnan(b)):
      b = np.copy(b_old)
      R = np.copy(R_old)
      break
    if cond_num < 1e-8 or np.isnan(cond_num):
      print 'Warning! Likelihood is flat.'
    # 2. Have we converged?
    if np.linalg.norm(b - b_old, np.inf) < epsilon:
      print 'Converged in %d steps.' % iter
      break
    
    # Update coefficients, covariance
    eta = offset + np.dot(X,b)
    mu = ilink_fun(eta)
    RI,Z1,Z2,Z3 = np.linalg.lstsq(R, np.eye(Np))
    C = np.dot(RI , RI.T)
    #C = np.multiply(C , s**2) # Right now: Assumes no over/under-dispersion (s=1)
    
  # Update stats structure
  st = stats()
  st.dfe = NT-Np
  #st.s = 1
  #st.estdisp = 0
  st.se = np.squeeze( np.sqrt( np.array(np.diag(C))) ) # Need se to be column vector
  st.t = np.divide(np.squeeze(np.array(b)), st.se)
  st.p = 2 * scipy.stats.norm.cdf( -np.abs(st.t) )
  st.weights = np.diag(Xw)
  return b,C,st

def rescale_isi(y, intensity):
    spike_ind = np.array([0]) # need 0 at beginning for integration
    spike_ind = np.append(spike_ind, np.nonzero(y)[0])
    num_isi = len(spike_ind)-1
    
    # Make sure we have enough spikes!
    if num_isi<3: raise ValueError('Too few data points!')

    # Integrate conditional intensity function
    z = np.zeros(num_isi)
    for j in range(num_isi):
      z[j] = np.sum(intensity[spike_ind[j]:spike_ind[j+1]])

    # Transform rescaled ISI, if desired
    # If RESCALING is 'exp', ISI should be uniform
    # If RESCALING is 'identity', ISI should be exponential
    if RESCALING=='exp': z = 1 - np.exp(-z)
    elif RESCALING=='identity': 0
    else:
      raise ValueError('Unknown rescaling %s' % RESCALING)
    return z

def kstest(y, intensity):
  z = rescale_isi(y, intensity)
  if RESCALING=='exp':
    def rs_cdf(x): return scipy.stats.uniform.cdf(x)
  else:
    def rs_cdf(x): return scipy.stats.expon.cdf(x, 1)
  ecdf = ECDF(np.sort(z))
  acdf = rs_cdf(ecdf.x)
  ks_stat = np.round( np.max(np.abs(acdf - ecdf.y)) , 3)
  ks_ci = np.round( 1.36/np.sqrt(len(z)), 3 )
  return ks_stat, ks_ci

def sim(m,d,p):
  print 'Simulation started...'
  N = 1
  T = d.len_time
  sim_dn = np.zeros([N,T])
  intensity = np.zeros([N,T])
  burnin = int(p.get_burnin())
  plotspline(p.covariates[1].knots, m.coeff[p.index[1]])
  lags1,yint = plotspline(p.covariates[1].knots, m.coeff[p.index[1]])
  lags2,yspa1 = plotspline(p.covariates[2].knots, m.coeff[p.index[2]])
  junk,yspa2 = plotspline(p.covariates[3].knots, m.coeff[p.index[3]])
  junk,yspa3 = plotspline(p.covariates[4].knots, m.coeff[p.index[4]])
  junk,yspa4 = plotspline(p.covariates[5].knots, m.coeff[p.index[5]])
  int_order = int(lags1[-1])
  spa_order = int(lags2[-1])
  for t in range(burnin,T):
    if np.mod(t,1e3)==0: print t*d.dt
    lambda_t = np.zeros(N)
    for n in range(N):
      #### VVV CAN MAKE THIS MORE FLEXIBLE
      lambda_t[n] = m.coeff[0]
      #####
      lambda_t[n] += np.sum(sim_dn[n, t-1-np.arange(0,int_order)] * yint)
      #lambda_t[n] += np.sum(np.multiply(d.dn[p.covariates[2].channel, t-1-np.arange(0,spa_order)] , yspa1))
      #lambda_t[n] += np.sum(np.multiply(d.dn[p.covariates[3].channel, t-1-np.arange(0,spa_order)] , yspa2))
      #lambda_t[n] += np.sum(np.multiply(d.dn[p.covariates[4].channel, t-1-np.arange(0,spa_order)] , yspa3))
      #lambda_t[n] += np.sum(np.multiply(d.dn[p.covariates[5].channel, t-1-np.arange(0,spa_order)] , yspa4))
    #for k in range(5):
    #  if n>k: lambda_t[n] += np.sum(sim_dn[neighbors(n)-k-1, t-1-np.arange(0,spa_order)] * yspa2)
    intensity[:,t] = np.exp(lambda_t)
    dn_t = (np.random.poisson.rvs(intensity[:,t])>0)
    sim_dn[:,t] = dn_t  
  print 'Simulation done!'
  return data(sim_dn, d.time),intensity

def image(C, cax, cmap='jet'):
  plt.figure()
  plt.pcolor(C, cmap=plt.get_cmap(cmap), vmin=cax[0], vmax=cax[1])
  plt.ylim([0,C.shape[0]])
  plt.xlim([0,C.shape[1]])
  plt.gca().invert_yaxis()