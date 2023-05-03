import numpy as np 

''' 
in this implementation the dimensionallity of enegy surface is arbitary
This implementation no longer works for 1D potential

'''

def parallel_temper(energy_function,x0,nrep,nsteps,dim,refreq,Temp,maxdr,**kwargs):
 
# initialize energies
# u=np.zeros((nrep))
# for rep in range(nrep):
#  u[rep] = energy_function(x[rep],**kwargs) # initial potential energy of the particle
# Follow the x-position of the particle in each replica
# set up an empty trajectory
 x = x0.copy()
 
 u = energy_function(x,**kwargs) # initial potential energy of the particle
 
 xTraj = np.zeros((dim,nrep,nsteps))
 acc1=np.zeros((nrep))
 acc2=np.zeros((nrep))
 rej=np.zeros((nrep))
 ex1c=0
 ex2c=0

 for step in range(nsteps):
    
        # Propose a replica exchange REFreq part of the time
        if (np.random.rand() < refreq):
            
          # pick a random replica
            rxch0=np.random.random_integers(0,nrep-1)   #
            if rxch0==0 : 
              rxch1=1               # if it is the lowest temp then only one nieghbor
            elif rxch0==nrep-1 : 
              rxch1=nrep-2     # if it is the highest temp then only neighbor
            else : 
             rxch1= rxch0+(np.random.random_integers(0,1)*2-1)   # it can be one of its niegbors 

            # Calculate beta differences and potential energy differences
            deltaBeta = 1.0/Temp[rxch1] - 1.0/Temp[rxch0]
            deltaU = u[rxch1]-u[rxch0]
            
            # Follow the Metropolis algorithm
            if (deltaBeta*deltaU > 0.0):
                print("exchange at step %i"%step)
                # Accept - thereforce, exchange the systems
                # Exchange x positions and energies  
                x[:,[rxch0,rxch1]]=x[:,[rxch1,rxch0]]    # this is syntactical numpy and will be handy when the number of reps is arbitary 
                u[[rxch0,rxch1]]=u[[rxch1,rxch0]]
                ex1c+=1
                
            # if deltaBeta*deltaU is positive, accept with the corresponding probability    
            elif (np.random.rand() < np.exp(deltaBeta*deltaU)):
                # Accept - thereforce, exchange the systems
                # Exchange x positions and energies  
                x[:,[rxch0,rxch1]]=x[:,[rxch1,rxch0]]    # this is syntactical numpy and will be handy when the number of reps is arbitary 
                u[[rxch0,rxch1]]=u[[rxch1,rxch0]]
                ex2c+=1
            
            # Otherwise - each replica remains with its own position and
            # potential energy            
        
        # Perform a regular displacement Monte Carlo move for each replica
        for rep in range(nrep):
            # Propose a trial move
            xTrial = x[:,rep] + maxdr[rep]*(np.random.rand(dim)-0.5)
            # Calculate the potential energy at the proposed position
            uTrial = energy_function(xTrial,**kwargs)
            # Calculate the difference in potential energy between current and
            # proposed positions
            deltaU = uTrial - u[rep]
            
            # Follow the Metropolis algorithm
            if (deltaU < 0.0):
                # If the proposed energy is lower - accept.
                x[:,rep] = xTrial
                u[rep] = uTrial
                acc1[rep]+=1
                
            elif (np.random.rand() < np.exp(-deltaU/Temp[rep])): 
                # If the proposed energy is higher - accept only with probablity
                # exp(-beta*deltaU)
                x[:,rep] = xTrial
                u[rep] = uTrial
                acc2[rep]+=1
            
            else:  # - stay at the same position, with the same potential energy.
                rej[rep]+=1
        
            # Save the position in the trajectory and the bath ID
            xTraj[:,rep,step] = x[:,rep]

        iter_to_screen=10
        alpha=1.0 # learning rate
        if ((step%int(np.round(nsteps/iter_to_screen)) == 0) and (step != 0) ) : 
           for rep in range(nrep):
            print("step: {:d} Temp: {:.3f} acc1: {:.3f} acc2: {:.3f} rej: {:.3f} maxdr: {:.3f} dev: {:.3f}".format(step, Temp[rep], 
                                                      acc1[rep]/step, acc2[rep]/step, rej[rep]/step, maxdr[rep], -alpha*maxdr[rep]*(rej[rep]/step - 0.5) ))

        if ((step%5 == 0) and (step != 0)) : 
           for rep in range(nrep):
          # deviation from desired acc should give shift magnitude - the learning rate  
          # below is similar to a finite difference method (eg. euler method)
          # we want to optimize maxdr given the constraint rej=0.5
              maxdr[rep] -=   ((rej[rep]/step) - 0.5)*alpha*abs(maxdr[rep]) 
          #    maxdr[rep] -=  0.02 
          #  else  :
          #    maxdr[rep] += 0.02
 
 for rep in range(nrep):
  print("step: {:d} Temp: {:.3f} acc1: {:.3f} acc2: {:.3f} rej: {:.3f}".format(step, Temp[rep], acc1[rep]/step, acc2[rep]/step, rej[rep]/step))
  print("exchange rates: ex1c: {:.3f}  ex2c: {:.3f} ".format(ex1c/step,ex2c/step))
  print("final stepsizes: ")
  print(maxdr)

 return xTraj


def bin_centers(bin_edges):
    return (bin_edges[1:]+bin_edges[:-1])/2.

def main():

 from energy_functions import  doubleWellPot as dw
 from matplotlib import pyplot as plt

# Set the initial state of all replicas
 nrep = 4                           # Number of replicas
 x=np.zeros((nrep))
 for rep in range(nrep):
  x[rep] = np.random.rand()-0.5              # initial x position of each replice

 
# plt.figure()
# plt.plot(X,Y)
# plt.show() 

 dim = 1   # this is a 1D potiental - the pt algoritham is modified to handle multidimensional PES
          # but for special case of 1D we need to modify the inputs and outputs for processing 

# Set simulation parameters
 nsteps = 1000 # Number of MC sweeps 10000
 refreq = 0.5  # 0.001       # 0.5 Fraction of the time to attempt replica exchange
 Temp = np.array([0.01,0.05,0.1,1.0])     # Temperature at each replica
 maxdr = np.array([0.5,1,5,10])     # Maximal displacement in each MC step, for each replica
 ID = np.arange(nrep)+1      # this is a bath ID used only to plot  how the baths exchange 

 xTraj = parallel_temper(dw,x,nrep,nsteps,dim,refreq,Temp,maxdr)

 xTraj=xTraj[0,:,:]  # only need tomake this move when dealing with a 1D potential

 hists=[]
 edges=[]
 boltzs=[]       
 for rep in range(nrep):
  x_hist, x_bin_edges = np.histogram(xTraj[rep],bins=50,density=False)  # these are the real the simualted distributions

 # Calculate the system's true boltzmann distribution
  boltz = np.exp( -dw(xTraj[rep])/Temp[rep])
  boltz = boltz/(np.sum(boltz)) #Normalize

  hists.append(x_hist)
  edges.append(x_bin_edges)
  boltzs.append(boltz) 
    
 X=np.linspace(-1.5,1.5,100)
 Y=dw(X)+0.2
 

 rgb=np.random.uniform(size=(nrep,3))
 fig, axes = plt.subplots(nrep,3, figsize=(12,8))
 for rep in range(nrep):
  axes[rep,0].plot(xTraj[rep],linestyle='-',color=(rgb[rep]))    # low temp bath 
  axes[rep,1].plot(bin_centers(edges[rep]),hists[rep]/float(nsteps),linestyle='-',color=(rgb[rep]),marker='.',alpha=0.6,markerfacecolor='none',markeredgewidth=1.5,markersize=5,label='x%i'%rep )
  axes[rep,1].plot(X,Y,'--k')
  axes[rep,2].plot(xTraj[rep],boltzs[rep],'.',color=(rgb[rep]))    # low temp bath 
 plt.show()


if __name__ == "__main__":
    main()
