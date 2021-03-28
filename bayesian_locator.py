import pymc3 as pm
import numpy as np
from scipy.stats import gaussian_kde

class BayesiantdoaPositionerMultiModal:
    """Class for carrying out Bayesian tdoa positioning.
    Requires at least 4 stations (to solve for x,y,v,t1).
    Able to incorporate two different wave speeds (multimodal positioning).
    """
    
    def __init__(self,
                 stations,
                 x_lim=1000,# maximum box size (m)
                 v_a_mu=350,# mean of velocity prior (m/s) (acoustic): 350
                 v_a_sd=50,# standard deviation of velocity prior (m/s) (acoustic): 50
                 t_a_sd=0.02,# standard deviation of observed values (s) (acoustic): 0.02
                 v_s_mu=400,# mean of velocity prior (m/s) (seismic): 400
                 v_s_sd=50,# standard deviation of velocity prior (m/s) (seismic): 50
                 t_s_sd=0.02):# standard deviation of observed values (s) (seismic): 0.02

        t_lim_a = np.sqrt(2)*x_lim/v_a_mu# resulting max tdoa_a value, used for t1 limit (s)
        t_lim_s = np.sqrt(2)*x_lim/v_s_mu# resulting max tdoa_a value, used for t1 limit (s)
        
        # check if well posed
        if len(stations)<4:
            print("WARNING: at least 4 stations usually required for bayesian tdoa positioning!")
        
        self.x_lim = x_lim
        self.stations = stations
        
        self.v_a_mu = v_a_mu
        self.v_a_sd = v_a_sd
        self.t_a_sd = t_a_sd
        self.t_lim_a = t_lim_a
        
        self.v_s_mu = v_s_mu
        self.v_s_sd = v_s_sd
        self.t_s_sd = t_s_sd
        self.t_lim_s = t_lim_s

        
    def sample(self, tdoa_a=None, tdoa_s=None, draws=2000, tune=2000, chains=4, init='jitter+adapt_diag', verbose=True):
        "Carry out Bayesian inference"
        
        x_lim = self.x_lim
        stations = self.stations
        
        v_a_mu = self.v_a_mu
        v_a_sd = self.v_a_sd
        t_a_sd = self.t_a_sd
        t_lim_a = self.t_lim_a
        
        v_s_mu = self.v_s_mu
        v_s_sd = self.v_s_sd
        t_s_sd = self.t_s_sd
        t_lim_s = self.t_lim_s
    
        
        
        if type(tdoa_a) != type(None):
            # assert correct number of observations
            if len(tdoa_a) != len(stations):
                raise Exception("ERROR: number of observations must match number of stations! (%i, %i)"%(len(tdoa_a), len(stations)))

        if type(tdoa_s) != type(None):
            # assert correct number of observations
            if len(tdoa_s) != len(stations):
                raise Exception("ERROR: number of observations must match number of stations! (%i, %i)"%(len(tdoa_s), len(stations)))

        with pm.Model() as model:
        
            # Priors
            x = pm.Uniform("x", lower=-600, upper=1200, shape=2)# prior on the source location (m)
            
            if type(tdoa_a) != type(None):
                v_a = pm.Uniform("v_a", lower=v_a_mu-v_a_sd, upper=v_a_mu+v_a_sd)
                t1_a = pm.Uniform("t1_a", lower=-0.5*t_lim_a, upper=t_lim_a)# prior on the time offset (s)
            
            if type(tdoa_s) != type(None):
                v_s = pm.Uniform("v_s", lower=v_s_mu-v_s_sd, upper=v_s_mu+v_s_sd)
                t1_s = pm.Uniform("t1_s", lower=-0.5*t_lim_s, upper=t_lim_s)# prior on the time offset (s)
            
            # Physics
            d = pm.math.sqrt(pm.math.sum((stations - x)**2, axis=1))# distance between source and receivers
            
            if type(tdoa_a) != type(None):
                t0_a = d/v_a# time of arrival (TOA) of each receiver
                t_a = t0_a-t1_a# time difference of arrival (tdoa_a) from the time offset
            
            if type(tdoa_s) != type(None):
                t0_s = d/v_s# time of arrival (TOA) of each receiver
                t_s = t0_s-t1_s# time difference of arrival (tdoa_s) from the time offset
            
            # Observations
            if type(tdoa_a) != type(None):
                Y_a_obs = pm.Normal('Y_a_obs', mu=t_a, sd=t_a_sd, observed=tdoa_a)# we assume Gaussian noise on the tdoa_a measurements
            if type(tdoa_s) != type(None):
                Y_s_obs = pm.Normal('Y_s_obs', mu=t_s, sd=t_s_sd, observed=tdoa_s)# we assume Gaussian noise on the tdoa_s measurements

                
            # Posterior sampling
            trace = pm.sample(draws=draws, tune=tune, chains=chains, target_accept=0.95, init=init)#, step=step)# i.e. tune for 1000 samples, then draw 5000 samples        
           
            summary = pm.summary(trace)
        
        mu = np.array(summary["mean"])
        sd = np.array(summary["sd"])
        
        if verbose:
            print("Percent divergent traces: %.2f %%"%(trace['diverging'].nonzero()[0].size / len(trace) * 100))
        
        return trace, summary, mu, sd
    
    def fit_xy_posterior(self, trace):
        """Helper function to estimate mu and sd of samples from a distribution,
        designed for when the tails of the distributions are large or non-zero"""
        print('calculating mean and standard deviation')
        mu = [np.mean(trace['x'][:,i]) for i in range(2)]
        sd = [np.std(trace['x'][:,i] - mu[i]) for i in range(2)]
        
        return mu, sd
    
    
    def forward(self, x, v=346):
        "predict time of flight for given source position"

        d = np.linalg.norm(self.stations-x, axis=1)
        t0 = d/v# time of flight values
        return t0