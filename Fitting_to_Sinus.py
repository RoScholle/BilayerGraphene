import Construct_H_0 as ConstructionSite
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import scipy.optimize as opt


# In[4]:
N_x = 20
N_y = 20
N_tot = N_x*N_y

def Set_Params(N_x_New = 20,N_y_New = 20):
    
    global N_x 
    global N_y
    global N_tot
    
    N_x = N_x_New
    N_y = N_y_New
    N_tot = N_x*N_y

plt.rcParams.update({"text.usetex": True,})

def R_x(theta):
    return np.matrix([[ 1, 0           , 0           ],
                   [ 0, np.cos(theta),-np.sin(theta)],
                   [ 0, np.sin(theta), np.cos(theta)]])
  
def R_y(theta):
    return np.matrix([[ np.cos(theta), 0, np.sin(theta)],
                   [ 0           , 1, 0           ],
                   [-np.sin(theta), 0, np.cos(theta)]])
  
def R_z(theta):
    return np.matrix([[ np.cos(theta), -np.sin(theta), 0 ],
                   [ np.sin(theta), np.cos(theta) , 0 ],
                   [ 0           , 0            , 1 ]])

    
    # In[17]:


def Cos_Fit(x,a,b,c,d):
    
    return c*np.cos(a*x+b)+d


def Schwebung(x,a,b,c,d,e):
    
    return c*np.cos(a*x+b)*np.cos(d*x+e)
    
    
def PlotMagn(GapVector,FitCos_x = [False,False,False],FitBeat_x = [False,False,False],FitCos_y = [False,False,False],FitBeat_y = [False,False,False],CosFrequencyx = [2*np.pi/10,2*np.pi/10,2*np.pi/10],CosFrequencyy= [2*np.pi/10,2*np.pi/10,2*np.pi/10],Beat_Frequenciesx = [[2*np.pi/80,2*np.pi/7.],[2*np.pi/80,2*np.pi/7.],[2*np.pi/80,2*np.pi/7.]],Beat_Frequenciesy = [[2*np.pi/80,2*np.pi/7.],[2*np.pi/80,2*np.pi/7.],[2*np.pi/80,2*np.pi/7.]],Fullparams = [],x_index = 0,y_index = 0,Index_to_fit = 3*N_x+1,QuickReturn = False):
    
    def Cos_Fit(x,a,b,c,d):
    
        return c*np.cos(a*x+b)+d


    def Schwebung(x,a,b,c,d,e):

        return c*np.cos(a*x+b)*np.cos(d*x+e)
    
    N_ups = GapVector[:N_tot]
    N_downs = GapVector[N_tot:2*N_tot]
    D_UpDowns = GapVector[2*N_tot:3*N_tot] + 1j*GapVector[3*N_tot:]
    
    Filling = N_ups+N_downs
    
    Phi_is = np.angle(D_UpDowns)
    M_is = np.sqrt((N_ups-N_downs)**2+4*np.abs(D_UpDowns)**2)
    Theta_is = np.arccos((N_ups-N_downs)/M_is)
    S_x = 2*np.real(D_UpDowns)
    S_y = 2*np.imag(D_UpDowns)
    S_z = N_ups-N_downs
    
    Coordinates = np.array([S_x,S_y,S_z])
    
    Rot_Mat = R_y(np.pi/2-Theta_is[0])@R_z(-Phi_is[0])
    Coordinates = Rot_Mat@Coordinates
    
    S_z,S_y,S_x = Coordinates
    S_x = np.array(S_x).reshape(N_tot)
    S_y = np.array(S_y).reshape(N_tot)
    S_z = np.array(S_z).reshape(N_tot)
    D_UpDowns = (S_x/2 + 1j*S_y/2.)
    N_ups = (Filling + S_z)/2
    N_downs = (Filling - S_z)/2
    
    Phi_is = np.angle(D_UpDowns)
    M_is = np.sqrt((N_ups-N_downs)**2+4*np.abs(D_UpDowns)**2)
    Theta_is = np.arccos((N_ups-N_downs)/M_is)
    
    
    Coordinates = np.array([S_x,S_y,S_z])
    
    Rot_Mat = R_z(-Phi_is[Index_to_fit]+np.pi/2)
    Coordinates = Rot_Mat@Coordinates
    
    S_z,S_y,S_x = Coordinates
    S_x = np.array(S_x).reshape(N_tot)
    S_y = np.array(S_y).reshape(N_tot)
    S_z = np.array(S_z).reshape(N_tot)
    #print(S_x[0],S_y[0],S_z[0])
    D_UpDowns = np.array(S_x/2 + 1j*S_y/2).reshape(N_tot)
    N_ups = np.array((Filling + S_z)/2).reshape(N_tot)
    N_downs = np.array((Filling - S_z)/2).reshape(N_tot)
    
    Phi_is = np.angle(D_UpDowns)
    M_is = np.sqrt((N_ups-N_downs)**2+4*np.abs(D_UpDowns)**2)
    Theta_is = np.arccos((N_ups-N_downs)/M_is)
    M_is = M_is.reshape(N_y,N_x)
    S_x = S_x.reshape(N_y,N_x)
    S_y = S_y.reshape(N_y,N_x)
    S_z = S_z.reshape(N_y,N_x)


# In[16]:


    MinusArray = np.ones((N_y,N_x))
    for i in range(N_x):
        for j in range(N_y):
            MinusArray[j,i] = (-1)**(i+j)



    if(QuickReturn == True):
        return MinusArray*S_x,MinusArray*S_y,MinusArray*S_z

    print("ALONG X-AXIS")
    
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(16, 4))
    
    data1 = (MinusArray*S_x)[y_index,:]
    data2 = (MinusArray*S_y)[y_index,:]
    data3 = (MinusArray*S_z)[y_index,:]
    
    xvals = np.array(range(len(data1)))
    xpoints = np.linspace(0,len(data1),500)
    P1 =[]

    ax[0].plot(xvals,data1)
    ax[1].plot(xvals,data2)
    ax[2].plot(xvals,data3)

    if(Fullparams == []):
        
        p01 = (Beat_Frequenciesx[0][0],0,0.2,Beat_Frequenciesx[0][1],0.)
        p02 = (Beat_Frequenciesx[1][0],0,0.2,Beat_Frequenciesx[1][1],0.)
        p03 = (Beat_Frequenciesx[2][0],0,0.2,Beat_Frequenciesx[2][1],0.)
        
    else:
        p01 = Fullparams[0]
        p02 = Fullparams[1]
        p03 = Fullparams[2]
    
    if(FitBeat_x[0]== True):
    
        params1,pcov = opt.curve_fit(Schwebung,xvals,data1,p0 = p01)
        ax[0].plot(xpoints,Schwebung(xpoints,params1[0],params1[1],params1[2],params1[3],params1[4]))
        P1.append(params1)
    if(FitBeat_x[1]== True):
        
        params2,pcov = opt.curve_fit(Schwebung,xvals,data2,p0 = p02)
        ax[1].plot(xpoints,Schwebung(xpoints,params2[0],params2[1],params2[2],params2[3],params2[4]))
        P1.append(params2)
    if(FitBeat_x[2]== True):
        
        params3,pcov = opt.curve_fit(Schwebung,xvals,data3,p0 = p03)
        ax[2].plot(xpoints,Schwebung(xpoints,params3[0],params3[1],params3[2],params3[3],params3[4]))
        P1.append(params3)


    if(Fullparams == []):
        
        p01 = (CosFrequencyx[0],0,0.05,0.)
        p02 = (CosFrequencyx[1],0,0.05,0.)
        p03 = (CosFrequencyx[2],0,0.05,0.)
    else:
        p01 = Fullparams[0]
        p02 = Fullparams[1]
        p03 = Fullparams[2]
        

    if(FitCos_x[0]== True):
    
        params1,pcov = opt.curve_fit(Cos_Fit,xvals,data1,p0 = p01)
        ax[0].plot(xpoints,Cos_Fit(xpoints,params1[0],params1[1],params1[2],params1[3]))
        P1.append(params1)
        
    if(FitCos_x[1]== True):
        
        params2,pcov = opt.curve_fit(Cos_Fit,xvals,data2,p0 = p02)
        ax[1].plot(xpoints,Cos_Fit(xpoints,params2[0],params2[1],params2[2],params2[3]))
        P1.append(params2)
    if(FitCos_x[2]== True):
        
        params3,pcov = opt.curve_fit(Cos_Fit,xvals,data3,p0 = p03)
        ax[2].plot(xpoints,Cos_Fit(xpoints,params3[0],params3[1],params3[2],params3[3]))
        P1.append(params3)

    
    

    
    ax[0].set_xlabel(r"$S_x$",fontsize = 14)
    ax[1].set_xlabel(r"$S_y$",fontsize = 14)
    ax[2].set_xlabel(r"$S_z$",fontsize = 14)
    
    
    plt.show()
    
    print("ALONG Y-AXIS:")
    
    
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(16, 4))
    
    data1 = (MinusArray*S_x)[:,x_index]
    data2 = (MinusArray*S_y)[:,x_index]
    data3 = (MinusArray*S_z)[:,x_index]
    
    xvals = np.array(range(len(data1)))
    
    ax[0].plot(xvals,data1)
    ax[1].plot(xvals,data2)
    ax[2].plot(xvals,data3)
    
    if(Fullparams == []):

        p01 = (Beat_Frequenciesy[0][0],0,0.2,Beat_Frequenciesy[0][1],0.)
        p02 = (Beat_Frequenciesy[1][0],0,0.2,Beat_Frequenciesy[1][1],0.)
        p03 = (Beat_Frequenciesy[2][0],0,0.2,Beat_Frequenciesy[2][1],0.)
    else:
        p01 = Fullparams[0]
        p02 = Fullparams[1]
        p03 = Fullparams[2]
    
    if(FitBeat_y[0]== True):
    
        params1,pcov = opt.curve_fit(Schwebung,xvals,data1,p0 = p01)
        ax[0].plot(xpoints,Schwebung(xpoints,params1[0],params1[1],params1[2],params1[3],params1[4]))
        P1.append(params1)
    if(FitBeat_y[1]== True):
        
        params2,pcov = opt.curve_fit(Schwebung,xvals,data2,p0 = p02)
        ax[1].plot(xpoints,Schwebung(xpoints,params2[0],params2[1],params2[2],params2[3],params2[4]))
        P1.append(params2)
    if(FitBeat_y[2]== True):
        
        params3,pcov = opt.curve_fit(Schwebung,xvals,data3,p0 = p03)
        ax[2].plot(xpoints,Schwebung(xpoints,params3[0],params3[1],params3[2],params3[3],params3[4]))
        P1.append(params3)


    if(Fullparams == []):
        
        p01 = (CosFrequencyx[0],0,0.05,0.)
        p02 = (CosFrequencyx[1],0,0.05,0.)
        p03 = (CosFrequencyx[2],0,0.05,0.)
        
    else:
        p01 = Fullparams[0]
        p02 = Fullparams[1]
        p03 = Fullparams[2]
        
    if(FitCos_y[0]== True):
    
        params1,pcov = opt.curve_fit(Cos_Fit,xvals,data1,p0 = p01)
        ax[0].plot(xpoints,Cos_Fit(xpoints,params1[0],params1[1],params1[2],params1[3]))
        P1.append(params1)

    if(FitCos_y[1]== True):
        
        params2,pcov = opt.curve_fit(Cos_Fit,xvals,data2,p0 = p02)
        ax[1].plot(xpoints,Cos_Fit(xpoints,params2[0],params2[1],params2[2],params2[3]))
        P1.append(params2)
    def Cos_Fit(x,a,b,c,d):
        
        return c*np.cos(a*x+b)+d
    
    
    def Schwebung(x,a,b,c,d,e):
        
        return c*np.cos(a*x+b)*np.cos(d*x+e)
    if(FitCos_y[2]== True):
        
        params3,pcov = opt.curve_fit(Cos_Fit,xvals,data3,p0 = p03)
        ax[2].plot(xpoints,Cos_Fit(xpoints,params3[0],params3[1],params3[2],params3[3]))
        P1.append(params3)



    
    ax[0].set_xlabel(r"$S_x$",fontsize = 14)
    ax[1].set_xlabel(r"$S_y$",fontsize = 14)
    ax[2].set_xlabel(r"$S_z$",fontsize = 14)
    
    
    plt.show()
    return MinusArray*S_x,MinusArray*S_y,MinusArray*S_z

# In[ ]:





# In[18]:
