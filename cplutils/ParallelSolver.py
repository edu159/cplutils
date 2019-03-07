import numpy as np
import itertools

class CouetteSolver:
    def __init__(self, phi0, dphi0, phi1, dphi1, dt, tini, tfin, dy, L, eta, rho, N=200):
        print tfin, tini, dt
        self.nsteps = int((tfin-tini)/dt) + 1
        self.tvec = np.linspace(tini, tfin, self.nsteps)
        self.ny = int(L/dy) + 1
        self.y = np.linspace(0.0, L, self.ny)
        if phi0 is not None:
            self.phi0 = np.array([phi0(t) for t in self.tvec])
            self.dphi0 = np.array([dphi0(t) for t in self.tvec])
        else:
            self.phi0 = np.zeros(self.nsteps)
            self.dphi0 = np.zeros(self.nsteps)

        if phi1 is not None:
            self.phi1 = np.array([phi1(t) for t in self.tvec])
            self.dphi1 = np.array([dphi1(t) for t in self.tvec])
        else:
            self.phi1 = np.zeros(self.nsteps)
            self.dphi1 = np.zeros(self.nsteps)

        self.dt = dt
        self.tfin = tfin
        self.tini = tini
        self.dy = dy
        self.L = L
        self.eta = eta
        self.rho = rho
        self.step = 0
        self.t = 0.0
        self.diff = eta/rho
        self.N = N
        self.u = np.zeros((self.nsteps, self.ny))
        self.stress = np.zeros((self.nsteps, self.ny))
        self.compute_stress = True
        self.compute_vel = True
        self.step_index = {}
        self.compute_all = False

    def update_phi0(self, ui):
        self.phi0[self.step] = ui
        if self.step > 0:
            self.dphi0  = np.gradient(self.phi0[0:self.step + 1], self.dt)

    def update_phi1(self, ui):
        self.phi1[self.step] = ui
        if self.step > 0:
            self.dphi1  = np.gradient(self.phi1[0:self.step + 1], self.dt)

    def _lamb(self, n):
        return n * np.pi / self.L

    def _Gu(self, t, n):
        tidx = int((t-self.tini)/self.dt)
        return 2*(self.dphi1[:tidx+1] * (-1.0)**n - self.dphi0[:tidx+1])/(np.pi*n)

    def _Fu(self, t, n):
        tidx = int((t-self.tini)/self.dt)
        return np.exp(-self.diff*(self._lamb(n))**2*(t-self.tvec[:tidx+1]))

    def compute_step(self, step=None):
        t = step * self.dt + self.tini
        s_u = s_stress = 0.0
        for n in range(1, self.N + 1):
            Cn = 2*(self.phi1[0] * (-1.0)**n - self.phi0[0])/(np.pi*n) # + integral(f(y))
            conv_int = np.trapz(self._Fu(t, n) * self._Gu(t, n), dx=self.dt)
            decay_term = np.exp(-self.diff*(self._lamb(n))**2*t)*Cn
            sum_terms = conv_int + decay_term
            if self.compute_vel:
                s_u += sum_terms*np.sin(self._lamb(n)*self.y)
            if self.compute_stress:
                s_stress += sum_terms*self._lamb(n)*np.cos(self._lamb(n)*self.y)
        if self.compute_stress:
            self.stress[step,:] = self.eta * ((self.phi1[step] - self.phi0[step])/self.L + s_stress)
        if self.compute_vel:
            self.u[step,:] = (self.phi1[step] - self.phi0[step])*(self.y/self.L) + self.phi0[step] + s_u


    def compute(self, times_in=[]):
        if not times_in:
            print "kj4", self.nsteps
            steps = range(self.nsteps)
            self.compute_all = True
        else:
            # TODO check if they are correct
            steps = list(((np.array(times_in) - self.tini)/self.dt).astype(int))
            steps.sort()
        for step in steps:
            print "{}%".format(float(step)/float(steps[-1])*100)
            self.compute_step(step)
        if self.compute_stress:
            self.stress = self.stress[steps,:]
        if self.compute_vel:
            self.u = self.u[steps,:]
        # Copy to not modify external times_in list
        if times_in:
            times = np.array(times_in)
        else:
            times = np.array(steps)*self.dt
        times.sort()
        self.step_index = {t : s for s,t in enumerate(times)}


    def cumm_stress(self, y0, y1, t=None):
        u0 = 0.0
        u1 = 0.0
        # yidx0 = int(y0/self.dy)
        # yidx1 = int(y1/self.dy)
        # if t is None:
        #     u0 = self.u[self.step, yidx0]
        #     u1 = self.u[self.step, yidx1]
        # else:
        #     u0 = self.compute_u(t)[yidx0]
        #     u1 = self.compute_u(t)[yidx1]
        # return  self.eta * (u1 - u0)
        #
    def next(self):
        self.step += 1
        self.t += self.dt




def lamb(n, L):
    return n * np.pi / L

def u_analy_0(y, t, N, L, eta, rho, A):
    v1 = A
    v0 = 0.0
    diff = eta/rho
    s_u = np.zeros(len(y))
    for n in range(1, N + 1):
        lambd = lamb(n, L)
        Cn = 2*(v1 * (-1.0)**n - v0)/(np.pi*n) # + integral(f(y))
        conv_int = 0.0
        decay_term = np.exp(-diff*(lambd)**2*t)*Cn
        sum_terms = conv_int + decay_term
        s_u += sum_terms*np.sin(lambd*y)
        # s_stress += sum_terms*lamb(n,L)*np.cos(lamb(n,L)*y)
    # stress = eta*(-A*(1-np.exp(-t/t0))/L + s_stress)
    u = A*y/L + s_u
    return u

def u_analy_1(y, t, N, L, eta, rho, A):
    # phi0 = A*t
    # phi1 = 0.0
    # f(x) = 0.0
    diff = eta/rho
    s = np.zeros(len(y))
    for n in range(1, N + 1):
        lambd = lamb(n, L)
        s += (np.exp(-diff * (lambd)**2*t) - 1) / n**3 * np.sin(lambd * y)
    return A*t*(1-y/L) + 2*A*L**2 / (np.pi**3*diff) * s

def u_analy_2(y, t, N, L, eta, rho, A, t0, v0, v1):
    # phi0 = A*(1-np.exp(-t/t0))
    # phi1 = 0.0
    # f(x) = 0.0
    diff = eta/rho
    s_u = np.zeros(len(y))
    s_stress = np.zeros(len(y))
    for n in range(1, N + 1):
        lambd = lamb(n, L)
        Cn = 2*(v1 * (-1.0)**n - v0)/(np.pi*n) # + integral(f(y))
        conv_int = 2*A/(np.pi*n*(1-diff*lambd**2*t0))*\
                   (np.exp(-t/t0) - np.exp(-diff*lambd**2*t))
        decay_term = np.exp(-diff*(lambd)**2*t)*Cn
        sum_terms = conv_int + decay_term
        s_u += sum_terms*np.sin(lambd*y)
        s_stress += sum_terms*lambd*np.cos(lambd*y)
    stress = eta*(-A*(1-np.exp(-t/t0))/L + s_stress)
    u = A*(1-np.exp(-t/t0))*(1-y/L) + s_u
    return u, stress

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    t0 = 160.0
    U = 2.0
    phi1= lambda t: 0.0
    phi0= lambda t: U*(1-np.exp(-t/t0))
    dphi0= lambda t: U/t0*np.exp(-t/t0)
    dphi1= lambda t: 0.0
    L =  20.6
    ncells = 100
    Nfourier = 100
    solver = CouetteSolver(phi0, dphi0,\
            phi1, dphi1,\
            dt=0.005, tini=0.0, tfin=400, dy=L/ncells, L=L,\
            eta=2.14, rho=0.81, N=Nfourier)
    u,stress = solver.compute([10000,30000, 50000,70000])# 5000, 80000])
    fig, (u_plot,stress_plot) = plt.subplots(1,2)
    # fig, u_plot = plt.subplots(1,1)
    u_plot.plot(solver.y/L, solver.u[10000,:]/U, 'r')
    u_plot.plot(solver.y/L, solver.u[30000,:]/U, 'g')
    u_plot.plot(solver.y/L, solver.u[50000,:]/U, 'b')
    u_plot.plot(solver.y/L, solver.u[70000,:]/U, 'm')
    # u_plot.plot(solver.y/L, solver.u[5000,:]/U, 'r')
    # u_plot.plot(solver.y/L, solver.u[40000,:]/U, 'g')
    # u_plot.plot(solver.y/L, solver.u[60000,:]/U,'b')
    # u_plot.plot(solver.y/L, solver.u[80000,:]/U,'m')
    u_plot.set_yticks(np.arange(0.0, 1.0, 0.1))
    # u_plot.set_xticks(np.arange(0.0, 1.0, 0.1))

    u_anal,s_anal = u_analy_2(solver.y, 10000*solver.dt, Nfourier,L, 2.14,0.81, U,t0, phi0(0), phi1(0))
    u_anal2,s_anal2 = u_analy_2(solver.y, 30000*solver.dt, Nfourier,L, 2.14,0.81, U,t0, phi0(0), phi1(0))
    u_anal3,s_anal3 = u_analy_2(solver.y, 50000*solver.dt, Nfourier,L, 2.14,0.81, U,t0,phi0(0), phi1(0))
    u_anal4,s_anal4 = u_analy_2(solver.y, 70000*solver.dt, Nfourier,L, 2.14,0.81, U,t0,phi0(0), phi1(0))
    u_plot.plot(solver.y/L, u_anal/U, 'ro')
    u_plot.plot(solver.y/L, u_anal2/U, 'go')
    u_plot.plot(solver.y/L, u_anal3/U, 'bo')
    u_plot.plot(solver.y/L, u_anal4/U, 'mo')

    stress_plot.plot(solver.y/L, -solver.stress[10000,:],'r')
    stress_plot.plot(solver.y/L, -solver.stress[30000,:],'g')
    stress_plot.plot(solver.y/L, -solver.stress[50000,:],'b')
    stress_plot.plot(solver.y/L, -solver.stress[70000,:],'m')
    stress_plot.plot(solver.y/L, -s_anal, 'ro')
    stress_plot.plot(solver.y/L, -s_anal2, 'go')
    stress_plot.plot(solver.y/L, -s_anal3, 'bo')
    stress_plot.plot(solver.y/L, -s_anal4, 'mo')
    stress_plot.set_yticks(np.arange(0.0, 0.2, 0.025))
    plt.show()
