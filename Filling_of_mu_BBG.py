#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import scipy.optimize as opt
import warnings
warnings.filterwarnings("ignore")


# In[2]:


###Hier systematischer:

def Band1(x,y):
    return np.real(-1/2 * ( ( 2 * ( tp )**( 2 ) + ( ( U )**( 2 ) + ( -2 * ( ( ( tp )**( 4 ) + ( ( t3 )**( 2 ) * ( 4 + ( t3 )**( 2 ) ) * ( np.abs( ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,-1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) ) )**( 4 ) + ( -4 * t3 * tp * ( ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,-1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) )**( 3 ) + ( ( -2 * ( -2 + ( t3 )**( 2 ) ) * ( tp )**( 2 ) + 4 * ( U )**( 2 ) ) * ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,-1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) * ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,-1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) + -4 * t3 * tp * ( ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,-1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) )**( 3 ) ) ) ) ) )**( 1/2 ) + 2 * ( 2 + ( t3 )**( 2 ) ) * (np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,-1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) * ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,-1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) ) ) ) )**( 1/2 ))

def Band2(x,y):
    return np.real(1/2 * ( ( 2 * ( tp )**( 2 ) + ( ( U )**( 2 ) + ( -2 * ( ( ( tp )**( 4 ) + ( ( t3 )**( 2 ) * ( 4 + ( t3 )**( 2 ) ) * ( np.abs( ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,-1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) ) )**( 4 ) + ( -4 * t3 * tp * ( ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,-1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) )**( 3 ) + ( ( -2 * ( -2 + ( t3 )**( 2 ) ) * ( tp )**( 2 ) + 4 * ( U )**( 2 ) ) * ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,-1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) * ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,-1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) + -4 * t3 * tp * ( ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,-1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) )**( 3 ) ) ) ) ) )**( 1/2 ) + 2 * ( 2 + ( t3 )**( 2 ) ) * ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,-1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) * ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,-1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) ) ) ) )**( 1/2 ))

def Band3(x,y):
    return np.real(-1/2 * ( ( 2 * ( tp )**( 2 ) + ( ( U )**( 2 ) + ( 2 * ( ( ( tp )**( 4 ) + ( ( t3 )**( 2 ) * ( 4 + ( t3 )**( 2 ) ) * ( np.abs( ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,-1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) ) )**( 4 ) + ( -4 * t3 * tp * ( ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,-1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) )**( 3 ) + ( ( -2 * ( -2 + ( t3 )**( 2 ) ) * ( tp )**( 2 ) + 4 * ( U )**( 2 ) ) * ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,-1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) * ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,-1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) + -4 * t3 * tp * ( ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,-1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) )**( 3 ) ) ) ) ) )**( 1/2 ) + 2 * ( 2 + ( t3 )**( 2 ) ) * ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,-1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) * ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,-1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) ) ) ) )**( 1/2 ))
def Band4(x,y):
    return np.real(1/2 * ( ( 2 * ( tp )**( 2 ) + ( ( U )**( 2 ) + ( 2 * ( ( ( tp )**( 4 ) + ( ( t3 )**( 2 ) * ( 4 + ( t3 )**( 2 ) ) * ( np.abs( ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,-1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) ) )**( 4 ) + ( -4 * t3 * tp * ( ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,-1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) )**( 3 ) + ( ( -2 * ( -2 + ( t3 )**( 2 ) ) * ( tp )**( 2 ) + 4 * ( U )**( 2 ) ) * ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,-1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) * ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,-1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) + -4 * t3 * tp * ( ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,-1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) )**( 3 ) ) ) ) ) )**( 1/2 ) + 2 * ( 2 + ( t3 )**( 2 ) ) * ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,-1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) * ( np.cos( ( 3 )**( -1/2 ) * y ) + ( 2 * np.cos( 1/2 * x ) * ( np.cos( 1/2 * ( 3 )**( -1/2 ) * y ) + complex( 0,-1 ) * np.sin( 1/2 * ( 3 )**( -1/2 ) * y ) ) + complex( 0,1 ) * np.sin( ( 3 )**( -1/2 ) * y ) ) ) ) ) ) )**( 1/2 ))


# In[3]:


M = 201
xs = np.linspace(-5,5,M)+0.0001
ys = np.linspace(-5,5,M)+0.0001

ll = (xs[0],ys[0])
ur = (xs[-1],ys[-1])

xs,ys = np.meshgrid(xs,ys)


# In[4]:


def honey(x,y):
    
    return Theta(-y+2*np.pi/np.sqrt(3))*Theta(y+2*np.pi/np.sqrt(3))*Theta(y+np.sqrt(3)*x +4*np.pi/np.sqrt(3))*Theta(y-np.sqrt(3)*x +4*np.pi/np.sqrt(3))*Theta(-y+np.sqrt(3)*x +4*np.pi/np.sqrt(3))*Theta(-y-np.sqrt(3)*x +4*np.pi/np.sqrt(3))

def Theta(x): 
    return np.heaviside(x,0.)


# In[5]:


def GetFilling(mu):
    
    Area1 = np.sum(honey(xs,ys)*Theta(Band1(xs,ys)+mu))*25/(M**2)
    Area2 = np.sum(honey(xs,ys)*Theta(Band2(xs,ys)+mu))*25/(M**2)
    Area3 = np.sum(honey(xs,ys)*Theta(Band3(xs,ys)+mu))*25/(M**2)
    Area4 = np.sum(honey(xs,ys)*Theta(Band4(xs,ys)+mu))*25/(M**2)
    
    return(Area1 + Area2 + Area3 + Area4)/ (2*AreaHex)


# In[6]:


AreaHex = np.sum(honey(xs,ys)*25/(M**2))


# In[8]:


P1 = np.array([-2/np.sqrt(3),2])*np.pi/np.sqrt(3)
P2 = np.array([2/np.sqrt(3),2])*np.pi/np.sqrt(3)
P3 = np.array([4/np.sqrt(3),0])*np.pi/np.sqrt(3)
P4 = np.array([2/np.sqrt(3),-2])*np.pi/np.sqrt(3)
P5 = np.array([-2/np.sqrt(3),-2])*np.pi/np.sqrt(3)
P6 = np.array([-4/np.sqrt(3),0])*np.pi/np.sqrt(3)

Points = np.array([P1,P2,P3,P4,P5,P6])
Points2 = [tuple(P1),tuple(P6),tuple(P5),tuple(P4),tuple(P3),tuple(P2),tuple(P1)]
Points_x = np.array([P1[0],P2[0],P3[0],P4[0],P5[0],P6[0]])
Points_y = np.array([P1[1],P2[1],P3[1],P4[1],P5[1],P6[1]])
plt.rcParams['contour.negative_linestyle'] = 'solid'


# In[11]:


def not_in_Hexagon(Vector):
    a,b = Vector[0],Vector[1]
    
    x,y = (B1*a + B2*b)/(np.linalg.norm(B1))

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


# In[12]:


def mask_outside_polygon(poly_verts, ax=None):
    """
    Plots a mask on the specified axis ("ax", defaults to plt.gca()) such that
    all areas outside of the polygon specified by "poly_verts" are masked.  

    "poly_verts" must be a list of tuples of the verticies in the polygon in
    counter-clockwise order.

    Returns the matplotlib.patches.PathPatch instance plotted on the figure.
    """
    import matplotlib.patches as mpatches
    import matplotlib.path as mpath

    if ax is None:
        ax = plt.gca()

    # Get current plot limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    # Verticies of the plot boundaries in clockwise order
    bound_verts = [(xlim[0], ylim[0]), (xlim[0], ylim[1]), 
                   (xlim[1], ylim[1]), (xlim[1], ylim[0]), 
                   (xlim[0], ylim[0])]

    # A series of codes (1 and 2) to tell matplotlib whether to draw a line or 
    # move the "pen" (So that there's no connecting line)
    bound_codes = [mpath.Path.MOVETO] + (len(bound_verts) - 1) * [mpath.Path.LINETO]
    poly_codes = [mpath.Path.MOVETO] + (len(poly_verts) - 1) * [mpath.Path.LINETO]

    # Plot the masking patch
    path = mpath.Path(bound_verts + poly_verts, bound_codes + poly_codes)
    patch = mpatches.PathPatch(path, facecolor='white', edgecolor='none',zorder = 2)
    patch = ax.add_patch(patch)

    # Reset the plot limits to their original extents
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    return patch


# In[13]:


N = 20

a1 = np.array([1.,0])
a2 = np.array([1/2.,np.sqrt(3)/2])

B1 = 2*np.pi*np.array([1,-1./np.sqrt(3)])
B2 = 4*np.pi/np.sqrt(3)*np.array([0,1.])

b1 = B1/N
b2 = B2/N

n_Array = np.linspace(0,N-1,N,dtype = int)
n,m = np.meshgrid(n_Array,n_Array)

Vxs = n*b1[0] + m*b2[0]
Vys = n*b1[1] + m*b2[1]


VecArray = np.stack([Vxs,Vys])
I = np.linalg.inv(np.column_stack((B1, B2)))

b1 = np.array([1,0])
b2 = np.array([0,1])
b3 = np.array([1,1])
b4 = np.array([0,-1])
b5 = np.array([-1,0])
b6 = np.array([-1,-1])


B_Array = np.array([b1,b2,b3,b4,b5,b6])

vec_new = np.zeros((2,N,N))
VecArrayNew = np.zeros((2,N,N))


c = 0
for i in range(N):
    for j in range(N):
        vec_new[:,i,j] = (I.dot(VecArray[:,i,j])+0.5)%1-0.5
        
        if(not_in_Hexagon(vec_new[:,i,j]) == True):
            c +=1
            for b in B_Array:
                if(not_in_Hexagon(vec_new[:,i,j]+b) == False):
                    vec_new[:,i,j] = vec_new[:,i,j]+b
                    break
                if(b[0] == b6[0] and b[1] == b6[1]):
                    print("STH WROOOONG")
                    
                    
        VecArrayNew[:,i,j] = B1*vec_new[0,i,j] + B2*vec_new[1,i,j]


# In[15]:
t0 = 3.16
tp = 0.381/t0
t3 = 0.38/t0
t4 = 0.#14/t0
U = 0.

def Plot_Fermi_Surface(U_New,Doping,tp_factor=4.,t3_factor = 4.):
   
    global t0 
    global tp
    global t3
    global t4
    global U
    
    t0 = 3.16
    tp = tp_factor*0.381/t0
    t3 = t3_factor*0.38/t0
    t4 = 0.#14/t0
    U = U_New
    
    n = 1.-Doping
    
    mu = opt.fsolve(lambda x: GetFilling(x) - n,-0.238,xtol=1.49012e-05)

    if(mu == -0.238):
        print("oops")

    fig,ax = plt.subplots(figsize = (8,8))
    ax.scatter(VecArrayNew[0,:,:],VecArrayNew[1,:,:])  

    ax.contour(xs,ys,Band1(xs,ys),levels = [mu],colors = ["red"])
    ax.contour(xs,ys,Band2(xs,ys),levels = [mu],colors = ["green"])
    ax.contour(xs,ys,Band3(xs,ys),levels = [mu],colors = ["blue"])
    ax.contour(xs,ys,Band4(xs,ys),levels = [mu],colors = ["purple"])


    mask_outside_polygon(Points2,ax)   

    for i in range(6):
        ax.plot([Points[i][0],Points[i-1][0]],[Points[i][1],Points[i-1][1]],"black",linestyle = "solid",linewidth = 2)

    ax.set_ylabel(r"$k_y$",fontsize = 25,rotation = "horizontal")
    ax.set_xlabel(r"$k_x$",fontsize = 25)
    ax.set_xlim(-4.5,4.5)
    ax.set_ylim(-4.5,4.5)

    ax.text(3,4,r"$D = {:.3}$".format(U),fontsize = 16)
    ax.text(3,3.5,r"$n = {:.3}$".format(GetFilling(mu)),fontsize = 16)
    ax.text(3,3.,r"$t_\perp = {:.2}$".format(tp),fontsize = 16)
    ax.text(3,2.5,r"$t_3 = {:.2}$".format(t3),fontsize = 16)


    plt.show()


# In[ ]:




