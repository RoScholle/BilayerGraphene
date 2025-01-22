import numpy as np
import scipy as sp
import scipy.optimize as opt
import matplotlib.pyplot as plt

n_filling = 0.
N_tot = 0
beta = 1.
U = 0.
T = 0.

def Construct_H_0(N_x_new,N_y_new,U_new=4.,t_prime = 0.,Doping = 0.,T_New=0.,Periodic_Boundaries_x = True,Periodic_Boundaries_y = True):
    
    global N_tot
    global n_filling
    global beta
    global U
    global H_0
    global T
    global N_x
    global N_y
    
    N_x = N_x_new
    N_y = N_y_new
    T = T_New
    U = U_new
    N_tot = N_x*N_y 
    
    V = 1.
    t = -V
    n_filling = 1-Doping#1 for half-filling 
    
    if(T == 0):
        beta = np.infty
    else:
        beta = 1./T
    
    H_0 = np.zeros([2*N_tot,2*N_tot])#Hamiltonian without any Deltas
    #GapVector = np.zeros(4*N_tot) #Contains first the n_ups, then n_downs, then Re(d_up d_down), then Im
    
    
    
    QuantumNumbers = []  #To write a single index insead of (x,y,s).
    for i in range(N_y):
        for j in range(N_x):
            for s in (1,-1):
                QuantumNumbers.append([j,i,s])
    
    
    
    
    H_0 = np.zeros([2*N_tot,2*N_tot],dtype = "complex")  #Fill up H_0 without the Gaps
    
    for i in range(2*N_tot):     
        for j in range(2*N_tot):  
            x1,y1,s1 = QuantumNumbers[i]
            x2,y2,s2 = QuantumNumbers[j]
            
            
            #Implement Nearest-Neighbour-Hopping:
            if(x1 == x2+1 or x1 == x2-1):
                if(s1 == s2 and y1 == y2):
                    H_0[i,j]+= t
            
            if(y1 ==y2+1 or y1 == y2-1):
                if(s1 == s2 and x1 == x2):
                    H_0[i,j]+= t
            
            
            #Implement Next-Nearest-Neighbour-Hopping:
            if(x1 == x2+1 or x1 == x2-1):
                if(y1 ==y2+1 or y1 == y2-1):
                    if(s1 == s2):
                    
                        H_0[i,j]+= t_prime
    
                    
                    
            if(Periodic_Boundaries_x == True): #implement boundaries
                
                #for nearest neighbour hopping
                if(x1 == 0 and x2 == N_x-1):
                    if(s1 == s2 and y1 == y2):
                        H_0[i,j]+= t
                if(x2 == 0 and x1 == N_x-1):
                    if(s1 == s2 and y1 == y2):
                        H_0[i,j]+= t
                        
            if(Periodic_Boundaries_y == True):
                
                #for nearest neighbour hopping
                if(y2 == 0 and y1 == N_y-1):
                    if(s1 == s2 and x1 == x2):
                        H_0[i,j]+= t
                if(y1 == 0 and y2 == N_y-1):
                    if(s1 == s2 and x1 == x2):
                        H_0[i,j]+= t
                        
                        
            if(Periodic_Boundaries_x == True):
                
                #For next nearest neighbour hopping
                if(x1 == 0 and x2 == N_x-1):
                    if(Periodic_Boundaries_y == True):
                        if(y1%N_y ==(y2+1)%N_y or y1%N_y == (y2-1)%N_y): #with modulus to account for corner-corner-hopping
                            if(s1 == s2):
                                H_0[i,j]+= t_prime
                    else:
                        if(y1 ==y2+1 or y1 == y2-1): #with modulus to account for corner-corner-hopping
                            if(s1 == s2):
                                H_0[i,j]+= t_prime
                                
                                
                if(x2 == 0 and x1 == N_x-1):
                    if(Periodic_Boundaries_y == True):
                        if(y1%N_y ==(y2+1)%N_y or y1%N_y == (y2-1)%N_y):
                            if(s1 == s2):
                                H_0[i,j]+= t_prime
                    else:
                        if(y1 ==y2+1 or y1 == y2-1):
                            if(s1 == s2):
                                H_0[i,j]+= t_prime
                            
            if(Periodic_Boundaries_y == True):
                
                #For next nearest neighbour hopping
                if(y2 == 0 and y1 == N_y-1):
                    if(x1 == x2+1 or x1 == x2-1):  # Here without the modulus, since else we would double count the corner-corner hopping
                        if(s1 == s2):
                            H_0[i,j]+= t_prime
                if(y1 == 0 and y2 == N_y-1):
                    if(x1 == x2+1 or x1 == x2-1):
                        if(s1 == s2):
                            H_0[i,j]+= t_prime
                    
            
    H_0 = np.array(H_0) 
    return H_0

def Get_H(GapVector,H_0): #Returns the full Hamiltonian for a given GapVector
    
    H = np.array(H_0)
    
    GapVector = np.real(GapVector)
    
    N_ups = GapVector[:N_tot]
    N_downs = GapVector[N_tot:2*N_tot]
    Re_up_down = GapVector[2*N_tot:3*N_tot]
    Im_up_down = GapVector[3*N_tot:4*N_tot]
    
    DiagonalEntries = U*np.append(N_downs,N_ups).reshape(2,N_tot).transpose().flatten()
    OffDiags = U*np.append(-Re_up_down+1j*Im_up_down,np.zeros(N_tot)).reshape(2,N_tot).transpose().flatten()[:-1]
    
    H = H + np.diag(DiagonalEntries) + np.diag(OffDiags,1) + np.diag(OffDiags.conjugate(),-1)
                  
    return np.array(H)

def n_occ(Energy,mu):
    return 1./(np.exp(beta*(Energy-mu))+1)



def Exp_Val_Matrix(Matrix,Energies,mu):
    N_Vector = n_occ(Energies,mu)

    BigM = (N_Vector*Matrix)@Matrix.transpose().conjugate()

    return BigM.transpose()
    


def Get_N_tot(mu,Energies):
    
    return np.sum(n_occ(Energies,mu))



def Next_Gap_Vector(GapVector,GetEnergy=False,ReturnBoth = False,ReturnSpectral = False,ReturnThree = False):

    H = Get_H(GapVector,H_0)
    epsilons, Eigenvectors = sp.linalg.eigh(H)
    
    
    N_filled = int(N_tot*n_filling)
    mu_initial = (epsilons[N_filled]+epsilons[N_filled-1])/2.
    
    mu = opt.root(lambda mu: Get_N_tot(mu,epsilons)-int(N_tot*n_filling),x0 = mu_initial).x[0]

    
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

        #To calculate the energy, Gunnarsson and Zaanen have a sign error in their very last term of eq. 2.
        #Check Energy in TD limit pls
        if(T == 0):
            Energy = sum(n_occ(epsilons,mu)*(epsilons)) - U*sum(N_ups*N_downs) + U*sum(np.real(D_UpDowns)**2)+U*sum(np.imag(D_UpDowns)**2)
        else:
            Energy = -T*np.sum(np.log(1+np.exp(-beta*(epsilons-mu)))) + mu*N_filled - U*sum(N_ups*N_downs) + U*sum(np.real(D_UpDowns)**2)+U*sum(np.imag(D_UpDowns)**2)
            
            
        if(ReturnSpectral ==True):
            return epsilons, mu
        if(ReturnBoth == True):
            return np.real(NewGapVector), Energy/int(N_tot)
        if(ReturnThree == True):
            return np.real(NewGapVector), Energy/int(N_tot), mu
        return Energy / int(N_tot) 
    return np.real(NewGapVector)


def Plot_Arrows(GapVector,N_x = 20,N_y = 20,enhance = 1.,Zoom = False,Print_Thetas = False,Show_x_z_plane = False):
    N_tot = N_x*N_y
    N_ups = GapVector[:N_tot]
    N_downs = GapVector[N_tot:2*N_tot]
    D_UpDowns = GapVector[2*N_tot:3*N_tot] + 1j*GapVector[3*N_tot:]


    Pol = N_ups-N_downs
    #Pol = Pol.reshape(N_y,N_x)
    Filling = N_ups+N_downs
    Filling = Filling.reshape(N_y,N_x)
    
    N_is = N_ups+N_downs
    Phi_is = np.angle(D_UpDowns)
    M_is = np.sqrt((N_ups-N_downs)**2+4*np.abs(D_UpDowns)**2)
    Theta_is = np.arccos((N_ups-N_downs)/M_is)
    M_is = M_is.reshape(N_y,N_x)
    
    if(Print_Thetas==True):
        print(Theta_is/np.pi)
    
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 4))

    if(Show_x_z_plane == True):
        for i in range(N_x):
            for j in range(N_y):

                Angle = Phi_is[N_x*j+i]-Phi_is[0] +np.pi/2
                dx = enhance *M_is.flatten()[N_x*j+i]*np.cos(Angle)*np.sin(Theta_is[N_x*j+i])
                dy = enhance *Pol.flatten()[N_x*j+i]
                axes[0].arrow(i-dx/2,j-dy/2,dx,dy,head_width = 0.12)
        
        
    else:
        for i in range(N_x):
            for j in range(N_y):

                Angle = Phi_is[N_x*j+i]-Phi_is[0] +np.pi/2
                dx = enhance *M_is.flatten()[N_x*j+i]*np.cos(Angle)
                dy = enhance *M_is.flatten()[N_x*j+i]*np.sin(Angle)
                axes[0].arrow(i-dx/2,j-dy/2,dx,dy,head_width = 0.12)
    axes[0].set_xlabel("Spin Amplitude and Orientation")
    if(Zoom == False):
        axes[0].set_xlim(-1,N_x)
        axes[0].set_ylim(-1,N_y)
    else:
        axes[0].set_xlim(-1,10.5)
        axes[0].set_ylim(-1,10.5)
    #print(Phi_is/np.pi)
    #print((Phi_is[1]-Phi_is[0])/(np.pi),(Phi_is[2]-Phi_is[1])/(np.pi),(Phi_is[3]-Phi_is[2])/(np.pi))

    im = axes[1].pcolormesh(Filling)
    axes[1].set_xlabel("Filling on each site")
    fig.tight_layout()
    #cax = fig.add_axes([0.99, 0.15, 0.02, 0.8])
    fig.colorbar(im,orientation='vertical')
    ###plt.savefig("N_18_20_U4_1over9Doping.pdf",format="pdf")
    
    return fig,axes
