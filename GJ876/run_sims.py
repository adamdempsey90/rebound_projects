#!/usr/bin/env python
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import rebound as rb
import migration_sim as mg






def pevolve(arg):
    tau_a = arg[0]
    tau_e = arg[1]
    d = arg[2]
    print '%d %d'%(tau_a,tau_e)
    params = {'tau_e': tau_e,
         'tau_a': tau_a,
         'abrupt': True,
          'integrator':'ias15',
          'direct':d,
         'dt': 0.01}

    planets=[]
    planets.append({'m': .37})
    planets.append({'m': 0.00082, 'rho': 1.04534352444679,'a':0.138,'e':0.01,'inc':0 })
    planets.append({'m':0.00256,'rho':1.04534352444679,'a': 0.25,'e': 0.01, 'inc':0 })
    planets.append({'m':5.0e-05,'rho':1.04534352444679,'a':0.45 ,'e':0.01,'inc':0 })
    fname='results_test/res_t%d_k%d_d%d'%(params['tau_e'],params['tau_a'],int(params['direct']))
    res,res1,sim = mg.evolve(planets,**params)
    res.save_state(fname=fname+'_a.dat')
    res1.save_state(fname=fname+'_b.dat')

    fig,axes=plt.subplots(7,2,figsize=(20,15))
    res.plot(fig=fig,axes=axes[:,0])
    res1.plot(fig=fig,axes=axes[:,1])
    axes[0,0].set_title('$\\tau_e = %d$, $\\tau_a=%d$'%(res.tau_e,res.tau_a),fontsize=15)
    fig.savefig(fname+'.png')
    return

ta_vals = np.array([3e1,1e2,3e2,1e3,3e3,1e4,3e4,1e5,3e5])
te_vals = np.array([3e1,1e2,3e2,1e3,3e3,1e4,3e4,1e5,3e5])
direct = np.array([True,False])
args=[]
for ta in ta_vals:
    for te in te_vals:
        for d in direct:
            args.append((ta,te,d))

pool = rb.InterruptiblePool()
pool.map(pevolve,args)
print 'Finished'


