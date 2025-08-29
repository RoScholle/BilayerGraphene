import numpy as np
import scipy as sp
import scipy.optimize as opt
import sys



N_x = int(sys.argv[1])           #Number of unit cells in x-direction of the underlying triangular lattice
N_y = int(sys.argv[2])            #Number of unit cells in y-direction of the underlying triangular lattice

t = 1.             #hopping strength. We recommend to set this to 1 and view it as an energy scale.
t_perp =  0.381/3.16      #nearest neighbour hopping strength.
t3 =  0.38/3.16
t4 = 1.*0.14/3.16

U = float(sys.argv[3])            #local Hubbard U (interaction strength). 3 is a moderate coupling

n_filling = float(sys.argv[4])     #The enforced average filling of the system. n_filling = 1 corresponds to half-filling.
                    #n_filling = 0.9 corresponds to 10% hole doping.

T =  float(sys.argv[5])           #Temperature T. We again recommend to set t=1, so T is in units of t.

Periodic_Boundaries = True #Apply periodic boundaries in x-direction
Displacement_Field =  float(sys.argv[6])




#From these, we calculate useful helping variables, i.e inverse temperature beta, number of sites N_tot and the Doping

N_tot = 4*N_x*N_y #Because there are 4 sites in a BBG unit cell
Doping = 1-n_filling
if(T==0):
    beta = np.infty
else:
    beta = 1./T



#We then construct the non-interacting(!) Hamiltonian H_0. In this code, this is done for the 
#non-interacting part of the Hubbard Hamiltonian on a hexagonal lattice with optional periodic boundaries 
#and an optional magnetic field.
#The code can be extended by replacing the matrix H_0 by a different matrix H_0 for a different model, e.g
#by adding additional next-to nearest neighbor hopping terms to the Hamiltonian. Note that t_prime only enters the
#equation in this non-interacting Hamiltonian H_0.
#Since it does not take much time, we construct the non-interacting Hamiltonian explicitely with for-loops
#in the hope of increasing readability.


H_0 = np.zeros([2*N_tot,2*N_tot],dtype = "complex") #This is the form of the Hamiltonian we intend to fill.
                                                    #It is a 2N_tot x 2N_tot Hamiltonian, since each site
                                                    #can have a spin up and a spin down electron.


QuantumNumbers = []     #We can see this as the labels of the rows and columns of the Hamiltonian.
                        #We can with this write a single index insead of (x,y,v,s).

for i in range(N_y):            #The order of these loops matters, because this way we enumerate the sites
    for j in range(N_x):        #from left to right an THEN from low to up. 
        for v in [0,1,2,3]:     #We label the 4 site in the BBG unit cell
            for s in (1,-1):    #We label up spins with 1 and down spins with -1.
                QuantumNumbers.append([j,i,v,s]) 


#Now we specify how to construct the non-interacting Hamiltonian. We loop though each entry H_0[i,j] of it:
for i in range(2*N_tot):     
    for j in range(2*N_tot):  
        x1,y1,v1,s1 = QuantumNumbers[i] #each index of H_0 corresponds to an x-coordinate, a y-coordinate,  
        x2,y2,v2,s2 = QuantumNumbers[j] #a site in the unit cell and a spin
        
        #Intra Layer Hoppings:
        
        ####Same Unit Cells:
        
        if(x1 == x2 and y1 == y2):
            if((v1 == 0 and v2 == 1) or (v1 == 1 and v2 == 0) or (v1 == 2 and v2 == 3) or (v1 == 3 and v2 == 2)):
                if(s1 == s2):
                    H_0[i,j]+= -t
       
    
        ####Same y Component:
        if(x2 == x1 + 1 and y2 == y1):
            if((v1 == 3 and v2 == 2) or (v1 == 1 and v2 == 0)):
                if(s1 == s2):
                    H_0[i,j]+= -t
                    
        if(x2 == x1 - 1 and y2 == y1):
            if((v1 == 2 and v2 == 3) or (v1 == 0 and v2 == 1)):
                if(s1 == s2):
                    H_0[i,j]+= -t
                
        
        ####Different y Component:
        if(x2 == x1 and y2 == y1 + 1):
            if((v1 == 1 and v2 == 0)):
                if(s1 == s2):
                    H_0[i,j]+= -t
                    
        if(x2 == x1 and y2 == y1 - 1):
            if((v1 == 0 and v2 == 1)):
                if(s1 == s2):
                    H_0[i,j]+= -t
        
        
        if(x2 == x1-1 and y2 == y1 + 1):
            if((v1 == 2 and v2 == 3)):
                if(s1 == s2):
                    H_0[i,j]+= -t
                                
        if(x2 == x1+1 and y2 == y1 - 1):
            if((v1 == 3 and v2 == 2)):
                if(s1 == s2):
                    H_0[i,j]+= -t                
                
        ####Inter Layer Hopping:  
        
        ####Nearest Neighbors:
        if(x2 == x1 and y2 == y1):
            if((v1 == 0 and v2 == 2) or (v1 == 2 and v2 == 0)):
                if(s1 == s2):
                    H_0[i,j]+= -t_perp                   
                
        
        ####t3
        if(Periodic_Boundaries == False):
            
            if(x2 == x1 and y2 == y1):
                if((v1 == 3 and v2 == 1) or (v1 == 1 and v2 == 3)):
                    if(s1 == s2):
                        H_0[i,j]+= -t3          

            if(x2 == x1 and y2 == (y1-1)):
                if((v1 == 3 and v2 == 1)):
                    if(s1 == s2):
                        H_0[i,j]+= -t3   

            if(x2 == (x1+1) and y2 == (y1-1)):
                if((v1 == 3 and v2 == 1)):
                    if(s1 == s2):
                        H_0[i,j]+= -t3   


            if(x2 == x1 and y2 == (y1+1)):
                if((v1 == 1 and v2 == 3)):
                    if(s1 == s2):
                        H_0[i,j]+= -t3   

            if(x2 == (x1-1) and y2 == (y1+1)):
                if((v1 == 1 and v2 == 3)):
                    if(s1 == s2):
                        H_0[i,j]+= -t3   
                
            
            #################Implement t4
            
            if(x2 == x1 and y2 == y1):
                if((v1 == 2 and v2 == 1) or (v1 == 1 and v2 == 2) or (v1 == 0 and v2 == 3) or (v1 == 3 and v2 == 0)):
                    if(s1 == s2):
                        H_0[i,j]+= -t4          

            if(x2 == x1-1 and y2 == y1):
                if((v1 == 2 and v2 == 1) or (v1 == 0 and v2 == 3)):
                    if(s1 == s2):
                        H_0[i,j]+= -t4   

            if(x2 == x1+1 and y2 == y1):
                if((v1 == 1 and v2 == 2) or (v2 == 0 and v1 == 3)):
                    if(s1 == s2):
                        H_0[i,j]+= -t4   

            if(x2 == x1-1 and y2 == y1+1):
                if((v1 == 0 and v2 == 3)):
                    if(s1 == s2):
                        H_0[i,j]+= -t4
                        
            if(x2 == x1+1 and y2 == y1-1):
                if((v1 == 3 and v2 == 0)):
                    if(s1 == s2):
                        H_0[i,j]+= -t4
                        
            if(x2 == x1 and y2 == (y1-1)):
                if((v1 == 2 and v2 == 1)):
                    if(s1 == s2):
                        H_0[i,j]+= -t4   

            if(x2 == x1 and y2 == (y1+1)):
                if((v1 == 1 and v2 == 2)):
                    if(s1 == s2):
                        H_0[i,j]+= -t4  
            
        if(Periodic_Boundaries == True): #implement boundaries
            
            #for intra layer nearest neighbour hopping
            #Periodic Boundaries in x direction:
            if(x1 == 0 and x2 == N_x-1):
                if(s1 == s2 and y1 == y2):
                    if((v1 == 0 and v2 == 1) or (v1 == 2 and v2 == 3)):
                        H_0[i,j]+= -t
                        
                        
            if(x1 == N_x-1 and x2 == 0):
                if(s1 == s2 and y1 == y2):
                    if((v1 == 1 and v2 == 0) or (v1 == 3 and v2 == 2)):
                        H_0[i,j]+= -t
                        
            if(x1 == 0 and x2 == N_x-1):
                if(s1 == s2 and y1 == y2-1):
                    if((v1 == 2 and v2 == 3)):
                        H_0[i,j]+= -t                        
            
            if(x1 == N_x-1 and x2 == 0 ):
                if(s1 == s2 and y1 == y2+1):
                    if((v1 == 3 and v2 == 2)):
                        H_0[i,j]+= -t 
            
            #Periodic Boundaries in y direction:
            if(y1 == 0 and y2 == N_y-1):
                if(s1 == s2 and x1 == x2):
                    if((v1 == 0 and v2 == 1)):
                        H_0[i,j]+= -t
                        
                        
            if(y1 == N_y-1 and y2 == 0):
                if(s1 == s2 and x1 == x2):
                    if((v1 == 1 and v2 == 0)):
                        H_0[i,j]+= -t
                        
                        
            if(y1 == 0 and y2 == N_y-1):
                if(s1 == s2 and x1 == x2-1):
                    if((v1 == 3 and v2 == 2)):
                        H_0[i,j]+= -t                        
            
            if(y1 == N_y-1 and y2 == 0):
                if(s1 == s2 and x1 == x2+1):
                    if((v1 == 2 and v2 == 3)):
                        H_0[i,j]+= -t   
                        
            #Correcting the corner-corner hopping:
            if(x1 == 0 and y1 == N_y-1 and x2 == N_x-1 and y2 == 0):
                if(s1 == s2):
                    if((v1 == 2 and v2 == 3)):
                        H_0[i,j]+= -t 
                        
            if(x2 == 0 and y2 == N_y-1 and x1 == N_x-1 and y1 == 0):
                if(s1 == s2):
                    if((v1 == 3 and v2 == 2)):
                        H_0[i,j]+= -t
                        
                        
            ##t_perp does not need PBC, so we continue with...
            
            ##t3:
            
            if(x2 == x1 and y2 == y1):
                if((v1 == 3 and v2 == 1) or (v1 == 1 and v2 == 3)):
                    if(s1 == s2):
                        H_0[i,j]+= -t3          

            if(x2 == x1 and y2 == (y1-1)%N_y):
                if((v1 == 3 and v2 == 1)):
                    if(s1 == s2):
                        H_0[i,j]+= -t3   

            if(x2 == (x1+1)%N_x and y2 == (y1-1)%N_y):
                if((v1 == 3 and v2 == 1)):
                    if(s1 == s2):
                        H_0[i,j]+= -t3   


            if(x2 == x1 and y2 == (y1+1)%N_y):
                if((v1 == 1 and v2 == 3)):
                    if(s1 == s2):
                        H_0[i,j]+= -t3   

            if(x2 == (x1-1)%N_x and y2 == (y1+1)%N_y):
                if((v1 == 1 and v2 == 3)):
                    if(s1 == s2):
                        H_0[i,j]+= -t3  
                        
            #######t4:
            
            if(x2 == x1 and y2 == y1):
                if((v1 == 2 and v2 == 1) or (v1 == 1 and v2 == 2) or (v1 == 0 and v2 == 3) or (v1 == 3 and v2 == 0)):
                    if(s1 == s2):
                        H_0[i,j]+= -t4          

            if(x2 == (x1-1)%N_x and y2 == y1):
                if((v1 == 2 and v2 == 1) or (v1 == 0 and v2 == 3)):
                    if(s1 == s2):
                        H_0[i,j]+= -t4   

            if(x2 == (x1+1)%N_x and y2 == y1):
                if((v1 == 1 and v2 == 2) or (v2 == 0 and v1 == 3)):
                    if(s1 == s2):
                        H_0[i,j]+= -t4   

            if(x2 == (x1-1)%N_x and y2 == (y1+1)%N_y):
                if((v1 == 0 and v2 == 3)):
                    if(s1 == s2):
                        H_0[i,j]+= -t4
                        
            if(x2 == (x1+1)%N_x and y2 == (y1-1)%N_y):
                if((v1 == 3 and v2 == 0)):
                    if(s1 == s2):
                        H_0[i,j]+= -t4
                        
            if(x2 == x1 and y2 == (y1-1)%N_y):
                if((v1 == 2 and v2 == 1)):
                    if(s1 == s2):
                        H_0[i,j]+= -t4   

            if(x2 == x1 and y2 == (y1+1)%N_y):
                if((v1 == 1 and v2 == 2)):
                    if(s1 == s2):
                        H_0[i,j]+= -t4 
##################################################################

        #####Displacement Fields

        if(x1 == x2 and y1 == y2):
            if(v1 == v2 == 0 or v1 == v2 == 1):
                if(s1 == s2):
                    H_0[i,j]+= Displacement_Field/2.
        
        if(x1 == x2 and y1 == y2):
            if(v1 == v2 == 2 or v1 == v2 == 3):
                if(s1 == s2):
                    H_0[i,j]-= Displacement_Field/2.   








#In the following, we define the functions that we will use to perform the self-consistency loop.
#Essential here is the concept of the GapVector: The GapVector will contain the 4N_tot self-consistent 
#fields, i.e. the Deltas. To keep the GapVector real valued, we decide to structure the gapvector like this:

#First, it contains the N_tot values of n_up, then the N_tot values of n_down, then the N_tot values of 
#the real part of the expectation value <c^dagger_{j,up} c_{j,down}> and then its imaginary parts.



def Get_H(GapVector,H_0): #Returns the full self-consistent interacting Hamiltonian for a given GapVector
    
    #"Read" out the GapVector 
    N_ups = GapVector[:N_tot] 
    N_downs = GapVector[N_tot:2*N_tot]
    Re_up_down = GapVector[2*N_tot:3*N_tot]
    Im_up_down = GapVector[3*N_tot:4*N_tot]
    
    #Insert the Deltas in their corresponding places in the Hamiltonian, multiplied by the interaction strength
    DiagonalEntries = U*np.append(N_downs,N_ups).reshape(2,N_tot).transpose().flatten()
    OffDiags = U*np.append(-Re_up_down+1j*Im_up_down,np.zeros(N_tot)).reshape(2,N_tot).transpose().flatten()[:-1]
    
    #putting everything together
    H = H_0 + np.diag(DiagonalEntries) + np.diag(OffDiags,1) + np.diag(OffDiags.conjugate(),-1)
                  
    return np.array(H)




#Here, we define the Fermi function shifted by the chemical potential, n_F(Energy-mu)
def n_occ(Energy,mu):
    return 1./(np.exp(beta*(Energy-mu))+1)



#Here, the "Matrix" argument will later be the Eigenvectors of the Hamiltonian.
#This function will return a matrix with all the expectation values <c^dagger_{j,sigma} c_{j',sigma'}>
#Works also at finite T, since we are using the occupation numbers.
def Exp_Val_Matrix(Matrix,Energies,mu):
    N_Vector = n_occ(Energies,mu)

    BigM = (N_Vector*Matrix)@Matrix.transpose().conjugate()

    return BigM.transpose()





#For a given set of single particle energies and a given chemical potential, this  function
#returns the average particle number. We need this function to determine the chemical potential
#in each iteration self-consistently.
def Get_N_tot(mu,Energies):
    
    return np.sum(n_occ(Energies,mu))




#One step of the iteration. This function takes a GapVector and returns the updated GapVector.
#Optionally, it also calculates the energy of state with the current GapVector and it can 
#return only the GapVector or also the energy or GapVector, Energy and Chemical Potential together.
#If one wants to return more than ust the GapVector, one must activate GetEnergy = True.
def Next_Gap_Vector(GapVector,GetEnergy=False,ReturnBoth = False,ReturnSpectral = False,ReturnThree = False):
    
    #First, get the full Hamiltonain
    H = Get_H(GapVector,H_0)
    #Diagonalize it to obtain all eigenvectors and eigenvalues
    #This is currently the clear bottleneck, scaling as N_tot^3!!! 
    #Note that ALL eigvals/eigvecs are needed, so Lanczos solvers often run into problems.
    epsilons, Eigenvectors = sp.linalg.eigh(H)
    
    
    #Here, we determine a chem. Pot. so that the filling is consistent.
    N_filled = N_tot*n_filling
    #Our initial guess for the chem. Pot. is the solution for zero T.
    mu_initial = (epsilons[int(N_filled+0.0001)]+epsilons[int(N_filled+0.0001)-1])/2.
    mu = opt.root(lambda mu: Get_N_tot(mu,epsilons)-N_filled,x0 = mu_initial).x[0]

    #We calculate the new expectation values and thus construct a new GapVector.
    N_ups_new = []
    N_downs_new = []
    Re_UpDown = []
    Im_UpDown = []
    ExpValMatrix = Exp_Val_Matrix(Eigenvectors,epsilons,mu)
    for i in range(N_tot):
        N_ups_new.append(np.real(ExpValMatrix[2*i,2*i]))
        N_downs_new.append(np.real(ExpValMatrix[2*i+1,2*i+1]))
        D_up_down = ExpValMatrix[2*i,2*i+1]
        Re_UpDown.append(np.real(D_up_down))
        Im_UpDown.append(np.imag(D_up_down))

    NewGapVector = np.append(N_ups_new,np.append(N_downs_new,np.append(Re_UpDown,Im_UpDown)))
    
    
    if(GetEnergy == True):
        N_ups = GapVector[:N_tot]
        N_downs = GapVector[N_tot:2*N_tot]
        D_UpDowns = GapVector[2*N_tot:3*N_tot] + 1j*GapVector[3*N_tot:4*N_tot]

        
        #This calculates the free energy of the system for T=0 and T!=0.
        #Note, that Gunnarsson and Zaanen have a sign error in their very last term of eq. 2.
        if(T == 0):
            Energy = sum(n_occ(epsilons,mu)*(epsilons)) - U*sum(N_ups*N_downs) + U*sum(np.real(D_UpDowns)**2)+U*sum(np.imag(D_UpDowns)**2)
        else:
            Energy = -T*np.sum(np.log(1+np.exp(-beta*(epsilons-mu)))) + mu*N_filled- U*sum(N_ups*N_downs) + U*sum(np.real(D_UpDowns)**2)+U*sum(np.imag(D_UpDowns)**2)
            
        if(ReturnSpectral == True):
            return epsilons, mu
        if(ReturnBoth == True):
            return np.real(NewGapVector), Energy/N_tot
        if(ReturnThree == True):
            return np.real(NewGapVector), Energy/N_tot, mu
        return Energy /(N_tot) 
    return np.real(NewGapVector)





#Here, we finally perform the iteration to let the system converge towards a final order. 
#We start by initializing the GapVector to a completely random configuration with roughly, 
#but not necessarily exactly, the right amount of electrons.


GapVector = np.append(n_filling*np.random.random(2*N_tot),np.random.random(2*N_tot)-0.5)

#To keep track of the Energies and the chemical Potential in each iteration, we write them
#into a list. It can be interesting to e.g. plot the Energies as a function of iteraton to
#see that they indeed decrease.
mus = []
Energies = np.array([])


#Here we iterate. We stop the iteration process in this example either after 1000 iterations
#or after none of the Deltas has changed by more than 10^-8. 
#For higher quality results, especially on bigger lattices, we recommend increasing the 
#maximum iteration number to e.g. 4000.

for i in range(1000):   
    

        
        
    NewGapVector,Energy,mu_new = Next_Gap_Vector(GapVector,GetEnergy = True,ReturnThree = True)
    Energies = np.append(Energies,Energy)
    mus.append(mu_new)
    if(max(np.abs(NewGapVector-GapVector))<10**(-8)):
        if(i>2):
            break
    Mixing = 0.4 + 0.1*np.random.random()   #We employ a mixing and add a small additional 
                                            #random mixing. Mixing significantly helps
                                            #the convergence process!
    GapVector = Mixing*GapVector + (1-Mixing)*NewGapVector



np.save("GapVector_Nx{}_Ny{}_U{}_Filling{:.8}_Temperature{}:.6_DispField{:.6}_BBG.npy".format(N_x,N_y,U,n_filling,T,Displacement_Field),NewGapVector)
np.save("Energies_Nx{}_Ny{}_U{}_Filling{:.8}_Temperature{:.6}_DispField{:.6}_BBG.npy".format(N_x,N_y,U,n_filling,T,Displacement_Field),Energies[-1])
np.save("GapVector_Nx{}_Ny{}_U{}_Filling{:.8}_Temperature{}:.6_DispField{:.6}_BBG.npy".format(N_x,N_y,U,n_filling,T,Displacement_Field),mu_new)










