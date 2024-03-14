# Suite of functions to use with the wrapped Fortran code.
# Enables parallelization of the calculations.

import numpy as np
from tqdm import tqdm
from joblib import Parallel, delayed 
import sys
from physical_definitions import *
import os
import json

import myLib # Fortran 90 suite of light-matter interaction codes

# arg_input = int(sys.argv[1])
# n_jobs=arg_input
n_jobs = 5
 
spontCase = 's'

conversion = 0.1482e-24 * 1.113e-16 # conversion from a.u. to S.I.

C3_au_conversion = (electronCharge*a0)**2/(4*np.pi*e0)

alpha_GS = 430 * conversion # polarizability

def simulation_wrapper(C3, initialCond,P,T,w0,titf,s0,lambd,Gamma,absProj,delta,alpha):
    
    return myLib.create_simulation_loss(C3, initialCond, P,T,w0,titf,s0,lambd,Gamma,absProj,delta, alpha, spontCase, 1)


def generateInitialCond(P,T,w0,n_samples=int(1e3), nAtoms=int(1), alpha=alpha_GS,lambd=lambd_trap):
    u0 = P*alpha / (np.pi*c*e0*w0**2)
    zR = np.pi*w0**2/lambd

    T = 0.1 * u0/kB

    omega_perp = np.sqrt(4*u0/(m*w0**2)) # radial trap frequency, Hz
    omega_par = np.sqrt(2*u0/(m*zR**2)) # longitudinal trap frequency, Hz  
    
    dx_par = np.sqrt(kB*T/(m*omega_par**2))
    dx_perp = np.sqrt(kB*T/(m*omega_perp**2))
    dv = np.sqrt(kB*T/m)
      
    vz0 = np.random.normal(loc=0, scale=dv, size=(n_samples, nAtoms))
    vy0 = np.random.normal(loc=0, scale=dv, size=(n_samples, nAtoms))
    vx0 = np.random.normal(loc=0, scale=dv, size=(n_samples, nAtoms))  
      
    [x0, y0, z0] = np.array([np.random.normal(loc=0, scale=1*dx_perp, size=(n_samples, nAtoms)),
                             np.random.normal(loc=0, scale=1*dx_perp, size=(n_samples, nAtoms)),
                             np.random.normal(loc=0, scale=1*dx_par, size=(n_samples, nAtoms))])
    
    initialCond = np.array([x0, y0, z0, vx0, vy0, vz0], order='F')
    initialCond = np.moveaxis(initialCond, 0, 2)
    
    return initialCond


def runSimulation_loss(C3,P,T,w0,titf, s0=[0], lambd=[lambd583], Gamma=[Gamma583], absProj=[[1.0,0,0]], delta=[0], alpha=[alpha_GS], 
                       n_samples=int(1e3)):
    print(n_jobs)
    
    initialCond = generateInitialCond(P,T,w0, n_samples, nAtoms=2, alpha=alpha_GS, lambd=lambd_trap)

    results = Parallel(n_jobs=n_jobs)(delayed(simulation_wrapper)(C3, initialCond[i], P,T,w0,titf,s0,lambd,Gamma,absProj,delta, alpha)
                                                for i in tqdm(range(n_samples)))
    
    # results = [simulation_wrapper(C3, nAtoms, initialCond[i,:,:], P, T, w0, titf, s0, lambd, Gamma, absProj, delta, alpha)
    #             for i in range(n_samples)]
    return results
        

def survivalProb(C3,P,T,w0,titf,s0=[0],lambd=[lambd583],Gamma=[Gamma583],absProj=[[1.0,0,0]],delta=[0],alpha=[alpha_GS],
                 n_samples=int(1e3)):
    # print(titf)
    results = runSimulation_loss(C3,P,T,w0,titf,s0,lambd,Gamma,absProj,delta,alpha,n_samples)
    #print('I am here')
    nAtoms = 2
    results = np.array(results)
    survProb = [sum(results[:,i])/len(results[:,i]) for i in range(nAtoms)]
    phScatt = list(results[:,-1]) # last index from results gives the photon scattering
    return [survProb, phScatt]

def convertPointToDash(value):
    valueName = str(value).split('.')
    return valueName

def cptd(value):
    return convertPointToDash(value)

if __name__ == "__main__":
    
    Gamma = Gamma583
    lambd = lambd583
    
    C3u = 0.26 * C3_au_conversion # converting from a.u. to S.I.
    C3g = 0.13 * C3_au_conversion # converting from a.u. to S.I.

    C3Vals = np.array([C3u, -C3g, 0]) * 1
    
    P = 0.6e-3 # trap power, in W
    T = 20e-6 # temperature, in K
    w0 = 0.5e-6 # trap waist, in um
    
    absProj = [[1.0, 0.0, 0.18]] # fixed absorption projection 
    alpha_E = 1*alpha_GS # valid since we have a magic condition for the transition in the lab
    
    t0 = 0

    # The following three variables should be always used with the given value (except in cases where they're being scanned)..
    # These are NOT the value used in the lab; rather, it's the value we should be using.
    # provided that we had a higher survival during imaging.
    
    tf = 1e-3
    dParam = 2.0
    s0 = 1.0
    
    dt = 1e-8
    titf = [t0, tf, dt]

    delta = -dParam*180e3 # detuning from resonance, in Hz
    
    n_samples = int(1e3) # we should do everything with the same # of samples. Choose either 3e2 or 1e3.
    
    wName = cptd(1e6*w0)
    pName = cptd(1e3*P)
    tfName = cptd(1e3*tf)
    dName = cptd(dParam)
    sName = cptd(s0)
    
    variableList = ["C3 (a.u.)", "Tweezer power (mW)", "w0 (um)", "s0", "Detuning (Gamma)", "exposure (ms)"] 
    
    #%% 
    ## The following functions are defined as different scans (i.e. "different measurements") ##
    
    def checkSurvP(P_list):
        fileName = "result_c3Fix_Pscan_w{0}-{1}_s{2}-{3}_d{4}-{5}_tf{6}-{7}".format(wName[0], wName[1],
                                                                                    sName[0], sName[1],
                                                                                    dName[0], dName[1],
                                                                                    tfName[0], tfName[1])
        variablesValues = [C3, list(P_list), w0, s0, dParam, tf]
        keyValues = zip(variableList, variablesValues)
        paramsDict = dict(keyValues)
        paramsDict['scanvar'] = 'Tweezer power (mW)'
        
        result = [survivalProb(C3, P_list[i],T,w0,titf,[s0],[lambd],[Gamma],absProj,[delta], [alpha_GS], n_samples)
                  for i in range(len(P_list))]
        return [result, fileName, paramsDict]

    def checkSurvW0_fixedU0(w0_list):
        fileName = "result_c3Fix_P{0}-{1}_wScan_s{2}-{3}_d{4}-{5}_tf{6}-{7}".format(pName[0], pName[1],
                                                                                    sName[0], sName[1],
                                                                                    dName[0], dName[1],
                                                                                    tfName[0], tfName[1])
        U0 = 150e-6 * kB
        P_fixedU0 = U0/alpha_GS * (np.pi * c * e0 * w0_list**2)
        variablesValues = [C3, P, list(w0_list), s0, dParam, tf]
        keyValues = zip(variableList, variablesValues)
        paramsDict = dict(keyValues)
        paramsDict['scanvar'] = 'w0 (um)'
        
        result = [survivalProb(C3, P_fixedU0[i],T,w0_list[i],titf,[s0],[lambd],[Gamma],absProj,[delta], [alpha_GS], n_samples)
                  for i in range(len(w0_list))]
        return [result, fileName, paramsDict]
		
    def checkSurvW0(w0_list):
        fileName = "result_c3Fix_P{0}-{1}_wScan_s{2}-{3}_d{4}-{5}_tf{6}-{7}".format(pName[0], pName[1],
                                                                                    sName[0], sName[1],
                                                                                    dName[0], dName[1],
                                                                                    tfName[0], tfName[1])
        variablesValues = [C3, P, list(w0_list), s0, dParam, tf]
        keyValues = zip(variableList, variablesValues)
        paramsDict = dict(keyValues)
        paramsDict['scanvar'] = 'w0 (um)'
        
        result = [survivalProb(C3,P,T,w0_list[i],titf,[s0],[lambd],[Gamma],absProj,[delta], [alpha_GS], n_samples)
                  for i in range(len(w0_list))]
        return [result, fileName, paramsDict]
        
    def checkSurvDeltaS0(s0_list, d_list):
        fileName = "result_c3Fix_P{0}-{1}_w{2}-{3}_sScan_dScan_tf{4}-{5}".format(pName[0], pName[1],
                                                                                 wName[0], wName[1],
                                                                                 tfName[0], tfName[1])
        variablesValues = [C3, P, w0, list(s0_list), list(d_list), tf]
        keyValues = zip(variableList, variablesValues)
        paramsDict = dict(keyValues)
        paramsDict['scanvar'] = ['s0', 'delta (Gamma)']
        delta_list = -d_list*180e3
        result = [[survivalProb(C3,P,T,w0,titf,[s0_list[i]],[lambd],[Gamma],absProj,[delta_list[j]], [alpha_GS], n_samples)
                  for i in range(len(s0_list))] for j in range(len(delta_list))] 
        return [result, fileName, paramsDict]
    
    
    def checkSurvC3(C3_list):
        fileName = "result_c3Scan_P{0}-{1}_w{2}-{3}_s{4}-{5}_d{6}-{7}_tf{8}-{9}".format(pName[0], pName[1],
                                                                                        wName[0], wName[1],
                                                                                        sName[0], sName[1],
                                                                                        dName[0], dName[1],
                                                                                        tfName[0], tfName[1])
        variablesValues = [list(C3_list/C3_au_conversion), P, w0, s0, dParam, tf]
        keyValues = zip(variableList, variablesValues)
        paramsDict = dict(keyValues)
        paramsDict['scanvar'] = 'C3 (a.u.)'
        
        result = [survivalProb(C3_list[i], P,T,w0,titf,[s0],[lambd],[Gamma],absProj,[delta], [alpha_GS], n_samples)
                  for i in range(len(C3_list))]
        return [result, fileName, paramsDict]
    
    def checkSurvDetuning(d_list):
        fileName = "result_c3_m0-13_P{0}-{1}_w{2}-{3}_s{4}-{5}_dScan_tf{6}-{7}".format(pName[0], pName[1],
                                                                                        wName[0], wName[1],
                                                                                        sName[0], sName[1],
                                                                                        tfName[0], tfName[1])
        variablesValues = [1, P, w0, s0, list(d_list), tf]
        keyValues = zip(variableList, variablesValues)
        paramsDict = dict(keyValues)
        paramsDict['scanvar'] = 'C3 (a.u.)'
        
        delta_list = d_list * 180e3
        
        result = [survivalProb(C3Vals, P,T,w0,titf,[s0],[lambd],[Gamma],absProj,[delta_list[i]], [alpha_GS], n_samples)
                  for i in tqdm(range(len(delta_list)))]
        return [result, fileName, paramsDict]
                
            
    Pscan = np.linspace(1, 201, 21)*1e-3 # scanning trap powers from 1 mW to 201 mW
    C3Scan = np.linspace(1,101,11)*C3_au_conversion # scanning C3 values from 1 to 10^3 x the value expected from Svetlana
    #C3Scan = C3Scan*0 
    w0Scan = np.linspace(0.4,2,7)*1e-6
    s0Scan = np.linspace(0.4,5,2)
    dScan = np.arange(-8,3+0.125/2,0.125)
    # dScan = np.array([-2])    
    # 
    #########################################################
    
    print(n_jobs)
    
    # 1D or 2D scan

    dim = 1

    if dim == 1:

        #result, fileName, paramsDict = checkSurvP(Pscan) # check survival for diff. P; all else is fixed
        #result, fileName, paramsDict = checkSurvC3(C3Scan) # check survival for diff. C3; all else is fixed
        # result, fileName, paramsDict = checkSurvW0_fixedU0(w0Scan) # check survival for diff. C3; all else is fixed
        result, fileName, paramsDict = checkSurvDetuning(dScan)
        
        surv1 = [result[i][0][0] for i in range(len(result))]
        surv2 = [result[i][0][1] for i in range(len(result))]
        phScatt = [result[i][1] for i in range(len(result))]
    
    elif dim == 2:
    
        result, fileName, paramsDict = checkSurvDeltaS0(s0Scan, dScan) # check survival for diff. C3; all else is fixed
        
        surv1 = [[result[i][j][0] for i in range(len(result))] for j in range(len(result[0]))]
        surv2 = [[result[i][j][1] for i in range(len(result))] for j in range(len(result[0]))]
        phScatt = [[result[i][j][2] for i in range(len(result))] for j in range(len(result[0]))]
        
   
    resultDict = {'dimension':dim, 'survival1': surv1, 'survival2': surv2, 'phScatt': phScatt}
    
    direct = 'C:\\Users\\c7041356\\Desktop\\LeonardoBG\\Science\\Project_TReqsSim\\code\\lightMatterMC_fortran90_v0102224\\results\\'
    
    os.chdir('results')
    
    # with open(fileName+'.json', "w") as savingFile:
    #     json.dump(resultDict, savingFile)
    
    # with open(fileName+'_params.json', "w") as savingParamFile:
    #     json.dump(paramsDict, savingParamFile)
    
    # os.chdir('..')
    
#os.remove(myLib.cpython-310-x86_64-linux-gnu.so) 
