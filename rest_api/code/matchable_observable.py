import numpy as np
import pandas as pd
import random as rand_

from math import pi, cos, sin, exp, floor, atan
from cmath import sqrt
from scipy.signal import StateSpace, impulse2, lsim2
from scipy.optimize import least_squares, dual_annealing, shgo
from scipy.linalg import lstsq


class MOLIwhiteBox:
    def __init__(self, method):
        self.C = np.array([[0., 0., 0., 1.]])
        if method is None:
            method = 'trivial'
        self.filter = method
        self.structure() #defines the model structure
        #self.describe() #print a model definition

    def structure(self):
        if self.filter == 'trivial':
            zeta = 0.
            omega = pi/12.0
            _p1 = -zeta*omega + omega*sqrt(zeta**2-1)
            _p2 = -zeta*omega - omega*sqrt(zeta**2-1)
            syscpoly = np.poly((-1.0*0.0353,_p1,_p2))
            sysc = np.concatenate((np.zeros([1]),np.flip(-syscpoly[1:])),0)

            A_ = np.concatenate((np.zeros([1,3]),np.identity(3)),0)
            A_ = np.concatenate((A_,sysc.reshape([4,1])),1)

            alphacand = np.poly(np.linalg.eig(A_)[0])
            A0 = np.concatenate((np.zeros([1,3]),np.identity(3)),0)
            A0 = np.concatenate((A0,-np.flip(alphacand[1:]).reshape([4,1])),1)
            self.A = A0

    def describe(self):
        print('This is a MOLI white Box alertness model')

    def fit(self, samples): #pipeline for simple fitting procedure
        self.awakeRegression(samples) #determine the parameters during the day
        self.parameters(samples) #determine the bio/physical parameters
        self.sim_Back_Forward(samples) #determine the homeostatic vicinities
        self.sleepRegression() #determine the night parameters with non linear least squares

    def awakeRegression(self, samples): #estimate the awake white box model parameters
        size_ = len(samples['levAl'])
        #create the model structure and the filtering state-space system
        A_, C_, I_ = np.transpose(self.A), np.transpose(self.C), np.identity(4)
        sFilt = StateSpace(A_, C_, I_, np.zeros([4,1]))
        #create the regressor information such as the B, output and regressor
        _B, regres, _y = np.zeros([4, size_]), np.array([]), []

        for i in range(size_):
            #create the alerteness filtered signal and the time data
            sigAl, _time = np.asarray(samples['levAl'][str(i)]), np.asarray(samples['time'][str(i)])

            #filter the input and output signals
            time_sim = _time - _time[0] #initialize time at t[0] = 0 (hours)
            impT, impResp = impulse2(sFilt, 0, time_sim, N = 10000) #simulate B info
            outT, outResp, x2 = lsim2(sFilt, sigAl, time_sim) #simulate L info

            #manipulate the regressor
            comp = np.zeros([impT.size,4*size_])
            comp[:,4*i:4*(i+1)] = impResp
            aux = np.concatenate((outResp,comp),1)

            if i == 0: #initialize the regressor, or add info to the regressor
                regres = aux
                _y = sigAl
            else:
                regres = np.concatenate((regres, aux), 0)
                _y = np.concatenate((_y, sigAl), 0)

        #theta = np.linalg.lstsq(regres, _y, rcond=None)[0] #solve the regression problem
        theta = lstsq(regres, _y, overwrite_a = True,
                      overwrite_b = True, lapack_driver = 'gelsy')[0] #'gelsy', 'gelss'

        #determine the parameters for each window
        for i in range(size_):
            init = 4*(i+1)
            final = 4*(i+1) + 4
            _B[:,i] = theta[init:final]
        self.L = theta[:4] #determine L
        self.B = _B #determine several B


    def parameters(self, samples): #determine the real biologycal parameters
        Ao = self.A + np.matmul(self.L.reshape([4,1]),self.C)
        self.w = (-Ao[2,3])**0.5 #determine the frequency
        self.tau = (self.w**2)/(-Ao[1,3]) #determine the discharge component

        #initialize the vectors
        size_ = len(samples['levAl'])
        DC, ho = np.zeros([size_, 1]), np.zeros([size_, 1])
        cphase, M = np.zeros([size_, 1]), np.zeros([size_, 1])

        #create the regressor matrix to retrive the B parameters
        _X = np.array([[self.w**2, 0., 1.0/self.tau],
                       [0., 1.0/self.tau, 1.0],
                       [1.0, 1.0, 0.]])

        for i in range(size_):
            DC[i] = self.B[0,i]*self.tau/(self.w**2) #determine the DC level
            #remove the DC information from the data
            compensator = np.array([-DC[i]*self.w**2, -DC[i]/self.tau, -DC[i]])
            _B = self.B[1:,i].reshape(3,1) + compensator
            #solve the linear equation problem
            _Th = np.linalg.lstsq(_X, _B, rcond=None)[0]

            ho[i] = _Th[0]
            tan_cphase = _Th[2]/(-self.w)/_Th[1]

            if tan_cphase < 0: #determine the phase considering possible algebraic
                cph = atan(tan_cphase) - pi
            else:
                cph = atan(tan_cphase)
            cphase[i] = cph - self.w*(samples['time'][str(i)][0] - 24*i)
            M[i] = _Th[1]/cos(cph)

        #estimate the parameters as the mean from each window
        self.M = np.mean(M)
        self.DC = np.mean(DC)
        self.ho = ho
        self.cphase = np.mean(cphase)

    def sim_Back_Forward(self, samples):
        #initialize the output data and the real A matrix from system
        self.hNight = pd.DataFrame(columns = {'h_s','h_w','t_s','t_w'})
        Ao = self.A + np.matmul(self.L.reshape([4,1]),self.C)
        _final, size_ = 0.0, len(samples['levAl'])

        for i in range(size_):
            _init = 24.0*i
            dt = _init - _final
            _final = _init + samples['windDecisions'][i]

            if i != (size_ - 1): #if the window is the last one, don't simulate forward
                step = (_final - samples['time'][str(i)][0]) / 100.0
                time_samp = [step*k + samples['time'][str(i)][0] for k in range(101)]
                time_samp = np.array(time_samp)
                #simulate the model from the first sampled moment to the sleeping moment
                sysF = StateSpace(Ao, self.B[:,i].reshape([4,1]), self.C, np.zeros([1,1]))
                outT, outRes = impulse2(sysF, T = (time_samp - time_samp[0]))
                circadian = self.M*cos(self.w*time_samp[-1] + self.cphase)
                #eliminate the circadian influence from the simulated signal
                self.hNight.loc[i,'h_s'] = outRes[-1] - circadian
                self.hNight.loc[i,'t_s'] = time_samp[-1]

            if i != 0: #if the window is the first one, don't simulate backwards
                step_n = (samples['time'][str(i)][0] - _init) / 100.0
                time_night = [step_n*k + _init for k in range(101)]
                time_night = np.array(time_night)
                #simulate the model from the first sampled moment to the awaking moment
                sysB = StateSpace(-Ao, self.B[:,i].reshape([4,1]), self.C, np.zeros([1,1]))
                #notice the negative A matrix --> backward simulation
                outT_n, outRes_n = impulse2(sysB, T = (time_night - time_night[0]))
                circadian = self.M*cos(self.w*time_night[0] + self.cphase)
                #eliminate the circadian influence from the simualted signal
                self.hNight.loc[i-1,'h_w'] = outRes_n[-1] - circadian
                self.hNight.loc[i-1,'t_w'] = time_night[0]

    def cost_function(self,x,sTime,sHom,wHom): #cost function for the non linear least squares
        c = [x[0]*(1-exp(-sTime[k]/x[1])) + sHom[k]*exp(-sTime[k]/x[1]) - wHom[k]
             for k in range(sTime.size)]
        return c

    def sleepRegression(self): #non linear least squares night parameters estimation
        #force the time to initialize at a zero reference
        sTime = np.asarray(self.hNight.loc[:,'t_w'].values - self.hNight.loc[:,'t_s'].values)
        #set the initial search state for the algorithm
        init__ = np.array([14.3, 2.6])
        aW, aS = self.hNight.loc[:,'h_w'].values, self.hNight.loc[:,'h_s'].values

        #estimate the parameters
        res_ = least_squares(self.cost_function, init__, args = (sTime, aS, aW) )
        self.sPar = res_
        self.y_ = res_.x[0]
        self.tau_e = res_.x[1]

    def predict(self, sTime, init_): #predict on especific moments provided a initial value

        alert, time_alert, decision = [], [], []
        for i in range(len(sTime['init'])):
            alSim, timeSim = np.array([]), np.array([])

            # determine the day time vector -- goes from the first sample until last
            step = (sTime['final'][i] - sTime['init'][i]) / 100.0
            _time_sim = np.array([sTime['init'][i] + k*step for k in range(101)])
            decision.append(_time_sim[-1])

            alSim = np.concatenate((alSim, self.lsim(_time_sim, init_)), 0)
            timeSim = np.concatenate((timeSim, _time_sim), 0)

            if i != (len(sTime['init'])-1):
                # determine the night time vector -- goes from last sample to first next day
                step = (sTime['init'][i+1] - sTime['final'][i]) / 100.0
                time_night = np.array([sTime['final'][i] + k*step for k in range(101)])
                # determine the circadian during night
                circ_night = np.asarray([self.M*cos(self.w*k + self.cphase) for k in time_night])
                #eliminate the influence of the circadian process
                h_o = alSim[-1] - circ_night[0]
                nightH = [self.y_*(1-exp(-dt/self.tau_e)) + h_o*exp(-dt/self.tau_e)
                         for dt in (time_night-time_night[0])]
                nightRes = np.asarray(nightH + circ_night)
                #reset the initial alertness level, to the last one from the night period
                init_ = nightRes[-1]
            else:
                nightRes = np.array([])
                time_night = np.array([])
            #concatenate the data for each window
            time_alert.append(list(np.concatenate((timeSim, time_night), 0)))
            alert.append(list(np.concatenate((alSim, nightRes), 0)))

        return {'decision': decision, 'levAl': alert, 'dt': time_alert}

    def lsim(self, Ts, init_): # simulate the dinamic system during the awake period
        #determine the B matrix
        w_t = Ts[0]*self.w + self.cphase
        k_1, k_2 = self.M*cos(w_t), -self.M*sin(w_t)*self.w
        h_o = init_ - k_1 - self.DC
        Bo = [[self.DC*self.w**2/self.tau],
              [k_2/self.tau + h_o*self.w**2 + self.DC*self.w**2],
              [(k_1 + self.DC + k_2*self.tau)/self.tau],
              [h_o + k_1 + self.DC]]
        #determine the A matrix
        Ao = self.A + np.matmul(self.L.reshape([4,1]), self.C)
        #simulate the system and return the result
        sysMOLI = StateSpace(Ao, Bo, self.C, np.zeros([1,1]))
        outT, outRes = impulse2(sysMOLI, T = (Ts - Ts[0]))
        return outRes

class forceEstimate:
    def __init__(self, alg):
        self.algorithm = alg
        self.omega = pi/12
        self.tau = 1/0.0353
        self.M = 2.52
        self.phi = -16.835*pi/12.0
        self.DC = 2.4
        self.y_ = 14.3
        self.tau_e = 1/0.381
        self.initStates = [self.omega,self.tau,self.M,self.phi,self.DC]

    def fit(self, samples):
        # Get the initial alertness value and initial parameters
        __init = samples['levAl']['0'][0]
        __initStates = self.initStates
        # Determine the lower and upper bounds
        lw = [pi/24, 18, 1.52, -2*pi, 1.4] # lower bounds
        up = [pi/6, 28, 3.52, 0, 3.4] # upper bounds
        # Estimate the parameters with a global optimization algorithm
        if self.algorithm == 'Annealing':
            res = dual_annealing(self.optSim, bounds = list(zip(lw, up)),
                            args = [samples, __init], x0 = __initStates)
        else:
            res = shgo(self.optSim, bounds = list(zip(lw,up)),
                            args = [samples, __init])
        # Determine the parameter vector theta
        theta = res.x
        self.omega, self.tau, self.M = theta[0], theta[1], theta[2]
        self.phi, self.DC =  theta[3], theta[4]
        # Estimate the night parameters
        awaking = self.simulate_Back_Forward(samples)
        self.sleepRegression(awaking)
        # Simulate the alertness model
        #return self.simulate(samples, __init)

    def optSim(self, par, samples, __init):
        # Simulate the system at known moments
        dataSim = self.sySim(par, samples, __init)
        # Get the comparison and real data for comparison
        compAl, realAl = dataSim['compAl'], dataSim['realAl']
        # Determine the error vector
        error = compAl.reshape([len(compAl), 1]) - realAl.reshape([len(realAl), 1])
        errorQuad = np.asarray([k*k for k in error]) # Quadratic error
        return np.sum(errorQuad) # Return the quadratic sumation error

    def sySim(self, par, samples, __init):
        size_ = len(samples['time'])
        # Initialize the parameters for simulation
        omega, tau, M, phi, DC = par[0], par[1], par[2], par[3], par[4]
        # Determine the initial homeostatic level and the initial time
        tZero = samples['time']['0'][0]
        initAl = __init - M*cos(omega * tZero + phi)
        # Simulate the system for each window
        for wii in range(size_):
            ind = str(wii)
            ##############   Simulate the day period   ########################
            # Detemine the especific times to evaluate
            timeGen = samples['time'][ind]
            timeRef = [timeSamp - tZero for timeSamp in timeGen]
            # Simulate the system at especific times
            wCirc = [M*cos(omega*k + phi) for k in timeGen]
            wHom = [(initAl - DC)*exp(-k/tau) + DC for k in timeRef]
            # Simulate the system forward until sleeping moment
            if samples['final'][wii] != timeGen[-1]:
                tSim = samples['final'][wii] - timeGen[0]
                __hInit = (initAl - DC)*exp(-tSim/tau) + DC
            else:
                __hInit = wHom[-1]

            #############   Simulate the night period   #######################
            if wii != size_ - 1: # Check if is different from last window
                period = samples['initial'][wii + 1] - samples['final'][wii]
                # Determine the homeostatic
                nightHom = self.y_*(1-exp(-period/self.tau_e)) + __hInit*exp(-period/self.tau_e)
                initAl = nightHom
                tZero = samples['initial'][wii + 1] # Set initial time
            # Check if is the first window to initialize the data set
            if wii != 0:
                aux = np.asarray(wCirc) + np.asarray(wHom)
                compAl = np.concatenate((compAl, aux), 0)
                realAl = np.concatenate((realAl,
                                        np.asarray(samples['levAl'][ind])), 0)
            else:
                compAl = np.asarray(wCirc) + np.asarray(wHom)
                realAl = np.asarray(samples['levAl'][ind])

        return {'compAl': compAl, 'realAl': realAl}

    def simulate_Back_Forward(self, samples):
        size_ = len(samples['levAl'])
        h = pd.DataFrame(columns = {'h_s','h_w','t_s','t_w'})
        # Simulate back towards the awaking moment
        for wii in range(size_):
            # Simulate backwards
            __init = samples['levAl'][str(wii)][0] # initial alertness
            __initTime = samples['time'][str(wii)][0] # initial alertness time
            initAl = __init - self.M*cos(self.omega * __initTime + self.phi)
            # Check if is different from the initial moment and if it is the
            # first window - first window do not simulate backward
            if (__initTime != samples['initial'][wii]) and (wii != 0):
                tSim = samples['initial'][wii] - __initTime
                h.loc[wii-1,'h_w'] = (initAl - self.DC)*exp(-tSim/self.tau) + self.DC
                h.loc[wii-1,'t_w'] = samples['initial'][wii]
            elif (__initTime == samples['initial'][wii]) and (wii != 0):
                h.loc[wii-1,'h_w'] = initAl
                h.loc[wii-1,'t_w'] = samples['initial'][wii]

            # Simulate forward
            __finalTime = samples['time'][str(wii)][-1] # final alertness time
            if (__finalTime != samples['final'][wii]) and (wii != size_ - 1):
                tSim = samples['final'][wii] - __initTime
                h.loc[wii,'h_s'] = (initAl - self.DC)*exp(-tSim/self.tau) + self.DC
                h.loc[wii,'t_s'] = samples['final'][wii]
            elif (__finalTime == samples['final'][wii]) and (wii != size_ - 1):
                __final = samples['levAl'][str(wii)][-1]
                __finalTime = samples['final'][wii]
                __circ = self.M*cos(self.omega * __finalTime + self.phi)
                h.loc[wii,'h_s'] = __final - __circ
                h.loc[wii,'t_s'] = __finalTime
        return h

    def sleepRegression(self, hNight):
        #non linear least squares night parameters estimation
        #force the time to initialize at a zero reference
        sTime = np.asarray(hNight.loc[:,'t_w'].values - hNight.loc[:,'t_s'].values)
        #set the initial search state for the algorithm
        init__ = np.array([14.3, 2.6])
        aW, aS = hNight.loc[:,'h_w'].values, hNight.loc[:,'h_s'].values
        #estimate the parameters
        res_ = least_squares(self.cost_function, init__, args = (sTime, aS, aW))
        self.sPar = res_
        self.y_ = res_.x[0]
        self.tau_e = res_.x[1]

    def cost_function(self,x,sTime,sHom,wHom):
        #cost function for the non linear least squares
        c = [x[0]*(1-exp(-sTime[k]/x[1])) + sHom[k]*exp(-sTime[k]/x[1]) - wHom[k]
             for k in range(sTime.size)]
        return c

    def simulate(self, samples, __init):
        size_ = len(samples['initial'])
        # Initialize the parameters for simulation
        omega, tau, M = self.omega, self.tau, self.M
        phi, DC = self.phi, self.DC
        # Determine the initial homeostatic level and the initial time
        tZero = samples['initial'][0]
        initAl = __init - M*cos(omega * tZero + phi)
        # Simulate the system for each window
        for wii in range(size_):
            ind = str(wii)
            ##############   Simulate the day period   ########################
            period = samples['final'][wii] - tZero
            # Detemine the especific times to evaluate
            timeRef = [k*period/100.0 for k in range(101)]
            timeGen = np.asarray(timeRef) + tZero
            # Simulate the system at especific times
            wCirc = [M*cos(omega*k + phi) for k in timeGen]
            wHom = [(initAl - DC)*exp(-k/tau) + DC for k in timeRef]
            #############   Simulate the night period   #######################
            # Determine the initial night homeostatic
            __hInit = wHom[-1]
            if wii != size_ - 1: # Check if is different from last window
                period = samples['initial'][wii + 1] - samples['final'][wii]
                # Create the zero referenced time vector - homeostatic
                nightRef = [k*period/100.0 for k in range(101)]
                # Create the non referenced time vector - circadian
                nightTime = np.asarray(nightRef) + samples['final'][wii]
                # Determine the circadian
                nightCirc = [M*cos(omega*k + phi) for k in nightTime]
                # Determine the homeostatic
                nightHom = [self.y_*(1-exp(-k/self.tau_e)) +
                            __hInit*exp(-k/self.tau_e) for k in nightRef]
                initAl = nightHom[-1]
                tZero = nightTime[-1] # Reinitialize the initial time
                nAl = np.asarray(nightCirc) + np.asarray(nightHom)
                nT = nightTime
            # Check if is the first window to initialize the data set
            dAl = np.asarray(wCirc) + np.asarray(wHom)
            dT = timeGen
            if (wii != 0) and (wii != size_ - 1):
                aux = np.concatenate((dAl, nAl), 0)
                Al = np.concatenate((Al, aux), 0)
                aux = np.concatenate((dT, nT), 0)
                time = np.concatenate((time, aux), 0)
            elif (wii == size_ - 1):
                Al = np.concatenate((Al, dAl), 0)
                time = np.concatenate((time, dT), 0)
            else:
                Al = np.concatenate((dAl, nAl), 0)
                time = np.concatenate((dT, nT), 0)

        return {'Al': Al.tolist(), 'time': time.tolist()}

class simAlert:
    def __init__(self, numDays):
        self.days = numDays
        self.time = pd.DataFrame()
        self.levAl = pd.DataFrame()
        self.levHom = pd.DataFrame()
        self.change = pd.DataFrame(columns = {'init','final'})
        self.omega = pi/12.0
        self.tau = 1/0.0353
        self.M = 2.52
        self.phi = -16.835*pi/12.0
        self.DC = 2.4
        self.y_ = 14.3
        self.tau_e = 1/0.381
        self.smpAl = pd.DataFrame(columns = {'time','levAl'} )

    def generate(self, windowDecision, resolution):
        initAl = 14.3
        self.resolution = resolution
        self.windows = windowDecision
        for i in range(self.days):
            unTime = [24*i + 24*j/self.resolution for j in range(resolution)]
            wTime = [j - 24*i for j in unTime if (j - 24*i) <= windowDecision[i]]
            sTime = [j - 24*i - windowDecision[i] for j in unTime
                     if (j - 24*i) >= windowDecision[i]]

            circadian = [self.M*cos(self.omega*k + self.phi) for k in unTime]
            wHom = [(initAl - self.DC)*exp(-k/self.tau) + self.DC for k in wTime]
            sHom = [self.y_*(1-exp(-k/self.tau_e)) + wHom[-1]*exp(-k/self.tau_e)
                    for k in sTime]
            Hom = np.concatenate((wHom,sHom))

            self.levHom.loc[:,i] = Hom[:resolution]
            self.levAl.loc[:,i] = circadian[:resolution] + Hom[:resolution]
            self.time.loc[:,i] = unTime[:resolution]
            initAl = Hom[-1]

    def randSample(self,ppH,seed):
        size_ = floor(max(self.windows)*ppH)
        time, levAl, day = [], [], []

        for i in range(self.days):
            step = floor(self.resolution/(24*ppH))
            final_ = floor(self.windows[i]*step*ppH)
            samples = range(0, final_, step)

            spH = floor(self.resolution/24.0)
            cruzeArr = [k - floor(spH/(2*ppH)) for k in range(0,floor(spH/ppH),1)]
            self.index = cruzeArr
            self.ind = samples
            cruzeInc = rand_.sample(cruzeArr, len(samples))
            if cruzeInc[0] < 0:
                cruzeInc[0] = 0
            if cruzeInc[-1] > 0:
                cruzeInc[-1] = 0
            samples_inc = np.asarray(cruzeInc) + samples

            time.append(list(self.time.loc[samples_inc,i].values))
            levAl.append(list(self.levAl.loc[samples_inc,i].values))

        self.smpAl.time = time
        self.smpAl.levAl = levAl
        samples = { 'time' : dict(self.smpAl.time),
                    'levAl' : dict(self.smpAl.levAl),
                    'windDecisions': self.windows}
        return (samples)
