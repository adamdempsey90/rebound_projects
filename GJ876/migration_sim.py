
import matplotlib.pyplot as plt
import numpy as np
import rebound as rb
import reboundx as rbx


def mod_2pi(ang):
    """ Mod an angle to be between 0 and 2*pi"""
    try:
        for i,p in enumerate(ang):
            while p < 0:
                p += 2*np.pi
            while p > 2*np.pi:
                p -= 2*np.pi
            ang[i] = p
    except TypeError:
        while ang < 0:
            ang += 2*np.pi
        while ang > 2*np.pi:
            ang-= 2*np.pi

    return ang

def mod_pi(ang):
    """ Mod an angle that's between -pi and pi to 0 to 2*pi"""
    try:
        for i,p in enumerate(ang):
            while p < -np.pi:
                p += 2*np.pi
            while p > np.pi:
                p -= 2*np.pi
            ang[i] = p
    except TypeError:
        while ang < -np.pi:
            ang += 2*np.pi
        while ang > np.pi:
            ang-= 2*np.pi

    return ang

def calc_timescale(t,tau,t0,ti):

    temp = ( 1- .5*(1 + np.tanh( (t - t0)/tau)) )/ti
    if temp < 1e-8:
        return np.infty
    else:
        return 1/temp


class Results():
    def __init__(self,times=np.linspace(0,100.,100),npart=3,fromfile=None,**init_pars):
        if fromfile is not None:
            self.read_state(fromfile)
        else:
            self.times = times
            self.npart = npart
            self.megno = np.zeros(times.shape)
            self.lyap = np.zeros(times.shape)
            self.energy = np.zeros(times.shape)
            self.avals = np.zeros((len(times),npart))
            self.evals = np.zeros((len(times),npart))
            self.ivals = np.zeros((len(times),npart))
            self.pvals = np.zeros((len(times),npart))
            self.lapvals = np.zeros(times.shape)
            self.dtvals = np.zeros(times.shape)
            self.ta_vals= np.inf*np.ones(times.shape)
            self.te_vals = np.inf*np.ones(times.shape)
            self.tend=times[-1]
            self.collision = False
            self.ejection = False
            self.restart_params = init_pars.copy()
            for key,val in init_pars.items():
                setattr(self,key,val)

            if self.integrator.lower() == 'whfast':
                self.integrator_i = 0
            else:
                self.integrator_i = 1


    def copy_final(self,res_old):
        self.avals[0,:] = res_old.avals[-1,:]
        self.evals[0,:] = res_old.evals[-1,:]
        self.ivals[0,:] = res_old.ivals[-1,:]
        self.pvals[0,:] = res_old.pvals[-1,:]
        self.lapvals[0] = res_old.lapvals[-1]
        self.megno[0] = res_old.megno[-1]
        self.lyap[0] = res_old.lyap[-1]
        self.energy[0] = res_old.energy[-1]
        self.dtvals[0] = res_old.dtvals[-1]
        self.ta_vals[0] = res_old.tau_a[-1]
        self.te_vals[0] = res_old.tau_e[-1]

        return
    def update(self,sim,i,chaos=False):

        if chaos:
            ps = sim.particles[:(self.npart+1)]
        else:
            ps = sim.particles
        self.avals[i,:] = [p.a for p in ps[1:]]
        self.evals[i,:] = [p.e for p in ps[1:]]
        self.ivals[i,:] = [p.inc for p in ps[1:]]
        self.pvals[i,:] = [p.P for p in ps[1:]]
        self.lapvals[i] = mod_2pi(ps[1].l - 3*ps[2].l + 2*ps[3].l)
        if chaos:
            self.megno[i] = sim.calculate_megno()
            self.lyap[i] = sim.calculate_lyapunov()

        self.energy[i] = sim.calculate_energy()
        self.dtvals[i] = sim.dt
        self.te_vals[i] = ps[3].tau_e if ps[3].tau_e is not np.nan else np.inf
        self.ta_vals[i] = ps[3].tau_a if ps[3].tau_a is not np.nan else np.inf

        return
    def cut(self):
        ind = self.times <= self.tend
        self.times = self.times[ind]
        self.avals = self.avals[ind,:]
        self.evals = self.evals[ind,:]
        self.pvals = self.pvals[ind,:]
        self.ivals = self.ivals[ind,:]
        self.lapvals = self.lapvals[ind]
        self.energy = self.energy[ind]
        self.dtvals = self.dtvals[ind]
        self.megno = self.megno[ind]
        self.lyap = self.lyap[ind]
        self.te_vals = self.te_vals[ind]
        self.ta_vals = self.ta_vals[ind]

        return
    def plot(self,fig=None,axes=None,fname=None,lap_lim=None,tstr=None):
        if axes is None:
            fig,axes=plt.subplots(7,1,sharex=True,figsize=(10,15))

        axa = axes[0]
        axe = axes[1]
        axl = axes[2]
        axE = axes[3]
        axmeg = axes[4]
        axper = axes[5]
        axmeg2 = axes[6]

        ind = self.times<=self.tend


        axa.plot(self.times[ind],self.avals[ind,0],'-b',
                 self.times[ind],self.avals[ind,1],'-g',
                 self.times[ind],self.avals[ind,2],'-r')
        axa.fill_between(self.times[ind],self.avals[ind,0]*(1-self.evals[ind,0]),
                self.avals[ind,0]*(1+self.evals[ind,0]),color='b',alpha=.2)
        axa.fill_between(self.times[ind],self.avals[ind,1]*(1-self.evals[ind,1]),
                self.avals[ind,1]*(1+self.evals[ind,1]),color='g',alpha=.2)
        axa.fill_between(self.times[ind],self.avals[ind,2]*(1-self.evals[ind,2]),
                self.avals[ind,2]*(1+self.evals[ind,2]),color='r',alpha=.2)
        axa.axhline(.349,color='r',linestyle='--')
        axa.axhline(.2186,color='g',linestyle='--')
        axa.axhline(.13598,color='b',linestyle='--')

        axe.plot(self.times[ind],self.evals[ind,0],'-b',
                 self.times[ind],self.evals[ind,1],'-g',
                 self.times[ind],self.evals[ind,2],'-r')

        axe.axhline(.075,color='r',linestyle='--')
        axe.axhline(.038,color='g',linestyle='--')
        axe.axhline(.252,color='b',linestyle='--')

        axl.plot(self.times[ind],mod_pi(self.lapvals[ind])*180./np.pi,'.')

        axE.plot(self.times[ind],abs((self.energy-self.energy[0])/self.energy[0])[ind])
        axmeg2.plot(self.times[ind],self.megno[ind])

        axmeg.plot(self.times[ind], ( 1/self.lyap )) # /(self.tau_e * self.K)) #/self.sim.tau_e)
        axper.plot(self.times[ind],self.pvals[ind,2]/self.pvals[ind,1],'-b',label='$P_3/P_2$')
        axper.plot(self.times[ind],self.pvals[ind,1]/self.pvals[ind,0],'-g', label='$P_2/P_1$')
        axa.set_ylabel('a',fontsize=15)
        axe.set_ylabel('e',fontsize=15)
        axl.set_ylabel('$|\\phi_{Laplace}|$',fontsize=15)
        axE.set_ylabel('$\\Delta E / E$',fontsize=15)
        axper.set_ylabel('$P_2/P_1$',fontsize=15)
        axmeg.set_ylabel('$\\tau_{lyap}(yrs)$',fontsize=15)
        axes[-1].set_xlabel('t (yrs)',fontsize=15)
        axper.legend(loc='upper right')
        axper.axhline(2,color='k')
        if lap_lim is not None:
            axl.set_ylim(-20,20)

        axmeg2.set_ylabel('MEGNO',fontsize=15)
        for ax in axes:
            ax.minorticks_on()

        if tstr is not None:
            axes[0].set_title(tstr,fontsize=15)
        if fname is not None:
            fig.savefig(fname)
        return

    def plot_laplace(self):
        ind2 = (self.times<=self.tend)
        t2 = self.times[ind2]

        vals = mod_pi(self.lapvals[ind2])*180./np.pi
        lp_mean =np.array([vals[t2<=x].mean() for x in t2])
        lp_std = np.array([vals[t2<=x].std() for x in t2])

        e3 = self.evals[ind2,2]
        e_mean =  np.array([ e3[t2<=x].mean() for x in t2 ])
        e_std = np.array([e3[t2<=x].std() for x in t2 ])

        e2 = self.evals[ind2,1]
        e2_mean =  np.array([ e2[t2<=x].mean() for x in t2 ])
        e2_std = np.array([e2[t2<=x].std() for x in t2 ])



        fig,axes=plt.subplots(4,1,figsize=(10,4),sharex=True)
        axes[0].plot(t2-t2[0],vals,'-b',alpha=.6)
        axes[0].plot(t2-t2[0],lp_mean,'-k',linewidth=3)
        axes[0].fill_between(t2-t2[0],lp_mean-lp_std,lp_mean+lp_std,color='r',alpha=.3)
        axes[1].plot(t2-t2[0],e3,'-b',alpha=.6)
        axes[1].plot(t2-t2[0],e_mean,'-k',linewidth=3)
        axes[1].fill_between(t2-t2[0],e_mean-e_std,e_mean+e_std,color='r',alpha=.3)
        axes[1].plot(t2-t2[0],e2,'-g',alpha=.6)
        axes[1].plot(t2-t2[0],e2_mean,'--k',linewidth=3)
        axes[1].fill_between(t2-t2[0],e2_mean-e2_std,e2_mean+e2_std,color='m',alpha=.6)

        axes[2].plot(t2-t2[0],1./self.lyap[ind2],'-b')
        axes[3].plot(t2-t2[0],self.megno[ind2],'-b')

        axes[0].set_ylabel('$\\phi$',fontsize=20)
        axes[1].set_ylabel('$e$',fontsize=20)
        axes[3].set_xlabel('$t - t_0$',fontsize=20)
        axes[2].set_ylabel('$\\tau_l$',fontsize=20)
        axes[3].set_ylabel('MEGNO',fontsize=15)
        plt.subplots_adjust(hspace=0)
        #axes[2].set_yscale('log')
        for ax in axes:
            ax.minorticks_on()
    def save_state(self,fname='results.dat'):

        with open(fname,'wb') as f:
            np.array([float(len(self.times)),float(self.npart)]).tofile(f)
            np.array([float(self.collision),float(self.ejection)]).tofile(f)
            np.array([self.tend,float(self.tau_e),float(self.tau_a),self.dt,float(self.integrator_i)]).tofile(f)
            self.times.tofile(f)
            self.megno.tofile(f)
            self.lyap.tofile(f)
            self.energy.tofile(f)
            self.avals.tofile(f)
            self.evals.tofile(f)
            self.ivals.tofile(f)
            self.pvals.tofile(f)
            self.lapvals.tofile(f)
            self.dtvals.tofile(f)
            self.te_vals.tofile(f)
            self.ta_vals.tofile(f)


        return
    def read_state(self,fname='results.dat'):
        dat = np.fromfile(fname)
        nt,npart,collision,ejection,self.tend,self.tau_e,self.tau_a,self.dt,self.integrator_i= dat[:9]
        npart = int(npart)
        nt = int(nt)
        self.collision = bool(collision)
        self.ejection = bool(ejection)
        self.npart = npart
        if self.integrator_i == 0:
            self.integrator = 'WHFast'
        else:
            self.integrator = 'ias15'
        dat = dat[9:]
        self.times = dat[:nt]
        dat = dat[nt:]
        self.megno = dat[:nt]
        dat = dat[nt:]
        self.lyap = dat[:nt]
        dat = dat[nt:]
        self.energy = dat[:nt]
        dat = dat[nt:]
        self.avals = dat[:nt*npart].reshape((nt,npart))
        dat = dat[nt*npart:]
        # avalsf
        self.ta_vals = dat[-nt:]
        dat = dat[:-nt]

        self.te_vals = dat[-nt:]
        dat = dat[:-nt]

        self.dtvals = dat[-nt:]
        dat = dat[:-nt]
        self.lapvals = dat[-nt:]
        dat = dat[:-nt]
        self.pvals = dat[-nt*npart:].reshape((nt,npart))
        dat = dat[:-nt*npart]
        self.ivals = dat[-nt*npart:].reshape((nt,npart))
        dat = dat[:-nt*npart]
        self.evals = dat[-nt*npart:].reshape((nt,npart))
        dat = dat[:-nt*npart]


        #self.evals = dat[:nt*npart].reshape((nt,npart))
        #dat = dat[nt*npart:]
        #self.ivals = dat[:nt*npart].reshape((nt,npart))
        #dat = dat[nt*npart:]
        #self.pvals = dat[:nt*npart].reshape((nt,npart))
        #dat = dat[nt*npart:]

        #self.lapvals = dat[:nt]
        #dat = dat[nt:]
        #self.dtvals = dat[:nt]
        #dat = dat[nt:]
        #self.te_vals = dat[:nt]
        #dat = dat[nt:]

        #self.ta_vals = dat[:nt]
        #dat = dat[nt:]


        return


def set_sim(init_params,dt=.02,integrator='WHFast'):
    sim = rb.Simulation()
    sim.units = ('yr','AU','Msun')
    for p in init_params:
            try:
                p.pop('rho')
            except KeyError:
                pass
            sim.add(**p)
    sim.move_to_com()
    sim.integrator = integrator
    ps = sim.particles
    if integrator == 'WHFast':
        sim.dt = 2*np.pi/ps[1].n * (1-ps[1].e)**(1.5) / 20
    else:
        sim.dt = dt
    sim.exit_max_distance = 5.
    sim.exit_min_distance = 3 * ps[2].a*np.sqrt(ps[2].m/3)


    return sim


def run_sim(sim,res,times,chaos=False,stopping_criterion=None):
    if chaos:
        sim.init_megno()

    for i,t in enumerate(times):
        try:
            sim.integrate(t,exact_finish_time=0)
        except rb.Escape as error:
            print error
            res.ejection = True
            res.tend = t
            res.cut()
            break
        except rb.Encounter as error:
            print error
            res.collision = True
            res.tend = t
            res.cut()
            break
        sim.calculate_orbits()
        res.update(sim,i,chaos=chaos)
        if stopping_criterion is not None:
            if stopping_criterion(sim,t):
                res.tend = t
                res.cut()
                print 'Stoping run at t=%f, stopping criterion met.'%res.tend
                return res
    res.tend = t
    if stopping_criterion is not None:
        print 'Finished before stopping criterion was met!'
    return res


def mmr_stop(sim,t):
# calculate_orbits() has been run
    ps = sim.particles
    outer_pair = ps[2].n/ps[3].n
    inner_pair = ps[1].n/ps[2].n
    if outer_pair <= 2.05:
        if inner_pair <= 2.0:
            return True
    return False
def a_stop(sim,t):
    ps = sim.particles
    if ps[3].a < 0.35:  # Outer planet
        return True
    return False

def run_constant_damping(sim,res,times,tau_e=1e2,tau_a=1e3,chaos=False,stopping_criterion=mmr_stop,direct=False):
    tau_a = res.tau_a
    tau_e = res.tau_e
    res.direct = res.direct

    rebx = rbx.Extras(sim)
    if direct:
        params = rebx.add_modify_orbits_direct()
    else:
        params = rebx.add_modify_orbits_forces()
    sim.particles[3].tau_a = -tau_a
    sim.particles[3].tau_e = -tau_e

    res = run_sim(sim,res,times,chaos=chaos,stopping_criterion=stopping_criterion)

    return res

def evolve(planets,**params):
    times = np.linspace(0,30*params['tau_a'],1e3)
    sim = set_sim(planets,integrator=params['integrator'],dt=params['dt'])
    res = Results(times,**params)
    res = run_constant_damping(sim,res,times,tau_e=params['tau_e'],tau_a=params['tau_a'],
            direct=params['direct'],stopping_criterion=mmr_stop,chaos=False)
    sim.particles[3].tau_a = np.infty
    sim.particles[3].tau_e= np.infty
    times1 = np.linspace(sim.t+sim.dt,sim.t+1e1*params['tau_a'],1e3)
    res1 = Results(times1,**params)
    res1 = run_sim(sim,res1,times1,chaos=True)

    return res,res1,sim

def para_evolve(args):
    planets = args[0]
    params = args[1]
    print 'dt = ',params['dt']
    res,res1,sim = evolve(planets,**params)
    return (res,res1)

def run_changing_damping(sim,res,times):
    ps = sim.particles
    tau_e_i = -ps[-1].tau_e
    K = abs(ps[-1].tau_a/tau_e_i)
    tot_time = times[-1]-times[0]
    tau = tot_time / 1e2
    t0 = times[0] + tot_time/2.
    tau_s = np.array([calc_timescale(t,tau,t0,tau_e_i) for t in times])
    plt.semilogy(times,tau_s,'-x')
    plt.ylabel('$\\tau_e$',fontsize=20)
    plt.xlabel('$t$',fontsize=20)
    for i,t in enumerate(times):
        #ps = sim.particles
        tau_e = calc_timescale(t,tau,t0,tau_e_i)
        ps[-1].tau_e = -tau_e
        ps[-1].tau_a = -tau_e * K
        try:
            sim.integrate(t*2*np.pi,exact_finish_time=0)
        except rb.Escape as error:
            print error
            res.tend = t
            break
        except rb.Encounter as error:
            print error
            res.tend = t
            break
        sim.calculate_orbits()
        res.update(sim,i,chaos=False)

    res.tend = t
    return res





