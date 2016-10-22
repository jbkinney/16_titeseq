import numpy as np
def hill(K,s,b,x):
#hill function
    return s*x/(K + x)+b


def make_x(K, amp, basal, x, fraction, sigmas, num_cells, boundaries):
    #makes x with 2 log-normal distributions. Kind of looks like real data
    y = hill(K,amp,basal,x)
    mu = np.log((y - (1 - fraction) * basal)/fraction) - sigmas[1]**2/2.
    mu0 = np.log(basal) - sigmas[0]**2/2.
    bound   = sigmas[1] * np.random.randn(int(np.round(num_cells*fraction)))+mu
    unbound = sigmas[0] * np.random.randn(int(np.round(num_cells*(1-fraction))))+mu0
    bound = np.hstack((bound, unbound))
    bound = np.exp(bound)
    x = np.histogram(bound, bins=boundaries)[0]
    return x, bound

def generate_bound(K, s, basal, fl, boundaries, sigmas, num_cells):
    #Generates a fake sort distribution for the data, truex, and then turns that into noisy sequencing data R
    truex = []
    R = []
    bounds = []
    for ii, f in enumerate(fl):
        x, bound = make_x(K, s, basal, f, 0.8, sigmas, num_cells, boundaries)
        bounds.append(bound)
        
    return bounds


def generate_counts(K, s, basal, k, S, T, fl, boundaries, sigmas):
    #Generates a fake sort distribution for the data, truex, and then turns that into noisy sequencing data R
    truex = []
    R = []
    bounds = []
    for ii, f in enumerate(fl):
        x, bound = make_x(K, s, basal, f, 0.8, sigmas, np.max([k*np.sum(S[ii]), 10]), boundaries)
        bounds.append(bound)
        x= np.array(x, dtype = float)/float(np.sum(x))
        x*=k/np.sum(x)
        truex.append(x)
        Sx = x * S[ii].sum()
        p = Sx / S[ii]
        p[p>1] = 1
        scale = 10**(np.random.rand()*4-2)
        #scale = 1
        R.append([float(np.random.binomial(int(T[ii][jj]/scale), p[jj])*scale) for jj in range(4)])
    
    truex = np.array(truex)
    R = np.array(R)
    return R, truex, bounds