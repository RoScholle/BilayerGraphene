import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
matplotlib.rcParams['text.usetex'] = True
import matplotlib.patches as mpat

# In[2]:

def Plotit(GapVector,N_x = 30,Quickreturn = False, returnaxes = False):
    
    def not_in_Hexagon(Vector):
        a,b = Vector[0],Vector[1]

        x,y = (b1*a + b2*b)/(np.linalg.norm(b1))

        P1 = np.array([-2/np.sqrt(3),2])/4
        P2 = np.array([2/np.sqrt(3),2])/4
        P3 = np.array([4/np.sqrt(3),0])/4
        P4 = np.array([2/np.sqrt(3),-2])/4
        P5 = np.array([-2/np.sqrt(3),-2])/4
        P6 = np.array([-4/np.sqrt(3),0])/4

        if(y>0.5):
            return True
        if(y<-0.5):
            return True
        if(y>(x-P6[0])*(P1[1]-P6[1])/(P1[0]-P6[0])):
            return True
        if(y<(x-P6[0])*(P5[1]-P6[1])/(P5[0]-P6[0])):
            return True
        if(y>(x-P3[0])*(P2[1]-P3[1])/(P2[0]-P3[0])):
            return True   


        if(y<(x-P3[0])*(P4[1]-P3[1])/(P4[0]-P3[0])):
            return True

        return False
    bb1 = np.array([1,0])
    bb2 = np.array([0,1])
    bb3 = np.array([1,1])
    bb4 = np.array([0,-1])
    bb5 = np.array([-1,0])
    bb6 = np.array([-1,-1])


    B_Array = np.array([bb1,bb2,bb3,bb4,bb5,bb6])

    N_y = N_x
    N_tot = N_x*N_y

    P1 = np.array([-2/np.sqrt(3),2])*np.pi/np.sqrt(3)
    P2 = np.array([2/np.sqrt(3),2])*np.pi/np.sqrt(3)
    P3 = np.array([4/np.sqrt(3),0])*np.pi/np.sqrt(3)
    P4 = np.array([2/np.sqrt(3),-2])*np.pi/np.sqrt(3)
    P5 = np.array([-2/np.sqrt(3),-2])*np.pi/np.sqrt(3)
    P6 = np.array([-4/np.sqrt(3),0])*np.pi/np.sqrt(3)

    Points = np.array([P1,P2,P3,P4,P5,P6])


    b1 = 2*np.pi*np.array([1,-1/np.sqrt(3)])
    b2 = 4*np.pi/np.sqrt(3)*np.array([0,1])
    l1 = b1/N_x
    l2 = b2/N_x

    a1 = np.array([1.,0])
    a2 = np.array([1./2.,np.sqrt(3)/2.])

    n_Array = np.linspace(0,N_x-1,N_x,dtype = int)
    n,m = np.meshgrid(n_Array,n_Array)

    Vxs = n*l1[0] + m*l2[0]
    Vys = n*l1[1] + m*l2[1]

    Cxs = n*a1[0] + m*a2[0]
    Cys = n*a1[1] + m*a2[1]

    VecArray = np.stack([Vxs,Vys])
    RealSpaceVecArray = np.stack([Cxs,Cys])
    I = np.linalg.inv(np.column_stack((b1, b2)))


    vec_new = np.zeros((2,N_x,N_x))
    VecArrayNew = np.zeros((2,N_x,N_x))
    c = 0
    for i in range(N_x):
        for j in range(N_x):
            vec_new[:,i,j] = (I.dot(VecArray[:,i,j])+0.5)%1-0.5

            if(not_in_Hexagon(vec_new[:,i,j]) == True):
                c +=1
                for b in B_Array:
                    if(not_in_Hexagon(vec_new[:,i,j]+b) == False):
                        vec_new[:,i,j] = vec_new[:,i,j]+b
                        break
                    if(b[0] == bb6[0] and b[1] == bb6[1]):
                        print("STH WROOOONG")


            VecArrayNew[:,i,j] = b1*vec_new[0,i,j] + b2*vec_new[1,i,j]
        
    def FourierTrafo(k,S):
    
        Component = np.einsum("...ij,ij",S,np.exp(1j*np.tensordot(k,RealSpaceVecArray,axes = 1)))    

        return Component

    def FullFourier(S):
        M = len(S)
        if(M == 1):
            Components = np.zeros((N_x,N_x),dtype = complex)
        else:
            Components = np.zeros((N_x,N_x,M),dtype = complex)
        for i in range(N_x):
            for j in range(N_x):
                k = VecArray[:,i,j]
                Component = FourierTrafo(k,S)
                Components[i,j] = Component

        return Components     
        
        

    N_ups = GapVector[:N_tot]
    N_downs = GapVector[N_tot:2*N_tot]
    D_UpDowns = GapVector[2*N_tot:3*N_tot] + 1j*GapVector[3*N_tot:]
    S_x = 2*np.real(D_UpDowns).reshape(N_x,N_x)
    S_y = 2*np.imag(D_UpDowns).reshape(N_x,N_x)
    S_z = (N_ups-N_downs).reshape(N_x,N_x)



    
    

    Components = FullFourier(np.array([S_x,S_y,S_z]))
    if(Quickreturn == True):
        return VecArrayNew,Components


    fig,ax = plt.subplots(figsize = (6,6))
    im = ax.scatter(VecArrayNew[0,:,:],VecArrayNew[1,:,:],c = np.linalg.norm(Components,axis = 2),marker = "H",s=8200*(np.sqrt(3)*2/N_x)**2,linewidth = 0,)  

    cax = fig.add_axes([0.94, 0.125, 0.05, 0.755])
    fig.colorbar(im, cax=cax, orientation='vertical')

    for i in range(6):
        ax.plot([Points[i][0],Points[i-1][0]],[Points[i][1],Points[i-1][1]],"black",linestyle = "dashed")
        ax.set_xlim(-4.5,4.5)
        ax.set_ylim(-4.5,4.5)
        ax.set_xlabel(r"$k_x$",fontsize = 20)
        ax.set_ylabel(r"$k_y$",fontsize = 20,rotation = "horizontal")
        
    if(returnaxes == True):
        return VecArrayNew,Components,fig,ax
        
    return VecArrayNew,Components

        


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

def Plot_Arrows(GapVector,N_x = 10,ShowSpin = True,ShowCharge = True,Scale = 0.5,SetMin = False,SetMax = False,Zoom = 1.):
    
    N_y = N_x
    N_tot = N_x*N_y
    
    
    N_ups = GapVector[:N_tot]
    N_downs = GapVector[N_tot:2*N_tot] 
    D_UpDowns = GapVector[2*N_tot:3*N_tot] + 1j*GapVector[3*N_tot:]


    


    
    if(ShowCharge == True):
        Filling = N_ups+N_downs
    
        if(SetMin == False):
            Minimum = min(Filling.flatten())
        else:
            Minimum = SetMin
        
        if(SetMax == False):
            Maximum = max(Filling.flatten())
        else:
            Maximum = SetMax
            
        Diff = Maximum - Minimum


        RelFil = (Filling.reshape(N_y,N_x)- Minimum)/Diff
        
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))

        Modulus = (N_x)*np.sqrt(3)  
        cmap = matplotlib.cm.get_cmap('plasma')




        for i in range(N_x):
            for j in range(N_y):

                xcoord = i*np.sqrt(3) + j*np.sqrt(3)/2 
                ycoord = 1.5*j

                color0 = cmap(RelFil[j,i])


                axes.add_patch(mpat.Circle(((xcoord+ 0.001)%Modulus - 0.001,ycoord),0.2,facecolor = color0))


        axes.set_xlim(-1,((N_x-1)*np.sqrt(3) + np.sqrt(3)/2 + 1)/Zoom)
        axes.set_ylim(-1,(N_y*1.5-0.5)/Zoom )
        
        axes.set_xlabel("Charges",fontsize = 16)
        
        
        newcax = inset_axes(axes,
                    width="7%",  # width = 50% of parent_bbox width
                    height="100%",  # height : 5%
                    bbox_to_anchor=(0.14, 0., 1, 1),
                    bbox_transform=axes.transAxes,
                    borderpad=0)
        
        plt.colorbar(matplotlib.cm.ScalarMappable(norm = matplotlib.colors.Normalize(Minimum,Maximum),cmap=cmap),cax = newcax)
        
        if(ShowSpin == False):
            plt.show()
        
    if(ShowSpin == True):
        
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(8,8))

        Modulus = (N_x)*np.sqrt(3)  

        Phi_is = np.angle(D_UpDowns)
        M_is = np.sqrt((N_ups-N_downs)**2+4*np.abs(D_UpDowns)**2)

        
        for i in range(N_x):
            for j in range(N_y):

                xcoord = i*np.sqrt(3) + j*np.sqrt(3)/2 
                ycoord = 1.5*j
                
                
                x1,y1 = (xcoord+ 0.001)%Modulus - 0.001, ycoord
                
                Angle = Phi_is[N_x*j+i]-Phi_is[0] +np.pi/2
                dx1 = 2*M_is[N_x*j+i]*np.cos(Angle)*Scale
                dy1 = 2*M_is[N_x*j+i]*np.sin(Angle)*Scale
                

                axes.arrow(x1-dx1/2,y1-dy1/2,dx1,dy1,head_width = 0.12)

        axes.set_xlim(-1,((N_x-1)*np.sqrt(3) + np.sqrt(3)/2 + 1)/Zoom)
        axes.set_ylim(-1,(N_y*1.5-0.5 )/Zoom)

        axes.set_xlabel("Spins",fontsize = 16)
        
        
        plt.show()
        
    return M_is



# In[ ]:




