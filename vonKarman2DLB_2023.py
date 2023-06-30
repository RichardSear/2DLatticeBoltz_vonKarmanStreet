'''
This is a modified version of a LB python code for 2D flow around obstacles by
# 
# Copyright (C) 2013 FlowKit Ltd, Lausanne, Switzerland
# E-mail contact: contact@flowkit.com
# 
# https://palabos.unige.ch/get-started/lattice-boltzmann/lattice-boltzmann-sample-codes-various-other-programming-languages/
# 
that was modified by Richard Sear summer 2020 to 2023
email: r.sear@surrey.ac.uk and richardsear.me

Below is the copyright text from the FlowKit code,
which also applies to this code:
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License, either
# version 3 of the License, or (at your option) any later version.
# 
'''
import numpy as np
import time
import matplotlib.pyplot as plt

'''
define some functions ....
 .... starting with LB Equilibrium distribution function.
'''
def equilibrium(rho,u):            
#    cu   = 3.0 * np.dot(c,u.transpose(1,0,2))
# not in Einstein summation, summing over i which is 0 to 8 and is over lattice vector,
# leaves array over j (x and y components) and over m and n which are lattice in x and y
# oddly setting optimize=True seems to slow code down a bit!
    lattice_vec_dot_u = 3.0 * np.einsum('ij,jmn',lattice_vec,u)#,optimize=True)
    usqr = 1.5*(u[0,:,:]**2+u[1,:,:]**2)
    feq = np.zeros((q,nx,ny))
    for i in range(0,q): 
        feq[i,:,:] = rho*t[i]*(1.0+lattice_vec_dot_u[i]+0.5*lattice_vec_dot_u[i]**2-usqr)
    return feq
###

"""
set up obstacle(s)
"""
def obstacle_def(r):
# cylinder centre at
    x0=nx/5
    y0=ny/2#+15
    print('cylinder at ',round(x0,2),round(y0,2))
##### now set True/False array for inside/outside disc
    inside_obs=np.zeros((nx,ny),dtype=bool)
    binary_obs=np.zeros((nx,ny))#,dtype=int)
    for i in range(0,nx):
        for j in range(0,ny):
            inside_obs[i,j]=False
            binary_obs[i,j]=0
# disc obstacle
            if( (i-x0)**2 + (j-y0)**2 < r**2  ): 
# inside obstacle
                inside_obs[i,j]=True
# eeeek!!!
#            if( (i-x0)**2 + (j-y0)**2 < (r+1.5)**2  ): 
#                binary_obs[i,j]=1
#            if( (i-x0)**2 + (j-y0+30)**2 < r**2  ): 
# inside obstacle
#                inside_obs[i,j]=True
                binary_obs[i,j]=1
#            if( i>x0 and i-x0<3.0*r and j==y0 ):
# inside 'fin' extension of obstacle
#                inside_obs[i,j]=True
#                binary_obs[i,j]=1               

#
    return inside_obs,binary_obs

'''
function to display while running or to save sequence of images
'''

def animation_plot(display_not_save_images):
    image1,image2=np.transpose(u[0,::n_coarse,::n_coarse])*(1.0-np.transpose(binary[::n_coarse,::n_coarse])), \
        np.transpose(u[1,::n_coarse,::n_coarse])*(1.0-np.transpose(binary[::n_coarse,::n_coarse]))
# update quiver plot
    im.set_UVC(image1,image2)
#
    vorticity = (np.roll(u[1,:,:], -1, axis=0) - np.roll(u[1,:,:], 1, axis=0)) - \
          (np.roll(u[0,:,:], -1, axis=1) - np.roll(u[0,:,:], 1, axis=1))
    print('min max vorticities ',np.amin(vorticity),np.amax(vorticity))
    im_vorticity.set_data(np.transpose(vorticity))
#
    stringy='Re ='+'{:.0f}'.format(round(Re, 0))+'   t/(2*r/u) = '+'{:.2f}'.format(round(tstep/T_d, 2))
    plt.title(stringy,fontsize=18)
#
    fig.canvas.draw()
    if(display_not_save_images):
       plt.pause(0.0001)
    else:
        plt.savefig("./pngs/ff"+str(int(1000000+tstep)).zfill(4)+".png")
#    
    return

# time run
tstart = time.time()

 
# lattice size, flow is along x axis, periodic boundary conditions along y
nx = 460
ny = 140
print('nx,  ny ',nx,ny)
r=20
print('radius of discs in lattice units ',round(r,2))
# works by setting desired Reynolds number, Re, then (below) setting kinematic viscosity
# to give this value of Re, at flow speed that is small, LB simulations are only stable
# at speeds << 1 in LB units. NB speed of sound is of order unity
Re = 100
print('Reynolds number for flow round cylinder (diameter as lengthscale) ',round(Re,2))
#
uBCLB     = 0.04                       # Velocity in lattice units.
print('x velocity imposed at LH boundary in LB units ',uBCLB,' must be << 1')
nuLB    = uBCLB*r/Re  #nulb=(2.0/omega-1.0)/6.0
print('kinematic viscosity nu in LB units ',round(nuLB,7))
omega = 1.0 / (3.0*nuLB+0.5); # Relaxation parameter.
tau=1.0/omega
print('LB relaxation rate omega = 1/tau ',round(omega,5))
print('relaxation time tau = ',round(tau,5))
# define reduced time units
T_d=2.0*r/uBCLB
print('reduced time diameter/imposed u = ',round(T_d,1),' LB steps')

###### Lattice Constants #######################################################
# standard D2Q9 LB
q = 9
#c = np.array([(x,y) for x in [0,-1,1] for y in [0,-1,1]]) # Lattice velocities.
#print(c)
#i=3
#jxy=1
#print(c[i,jxy])
lattice_vec=np.zeros((9,2),dtype=int)
lattice_vec[0]=[0,0]
lattice_vec[1]=[0,-1]
lattice_vec[2]=[0,1]
lattice_vec[3]=[-1,0]
lattice_vec[4]=[-1,-1]
lattice_vec[5]=[-1,1]
lattice_vec[6]=[1,0]
lattice_vec[7]=[1,-1]
lattice_vec[8]=[1,1]
for i in range(q):
    print('lattice velocity vector ',i,' = ',lattice_vec[i])
# set weights, off diagonal ones are 1/36, all for D2Q9 model
t = np.ones(q)/36.0
# 4 cardinal point weights 1/9
t[1]=1.0/9.0
t[2]=1.0/9.0
t[3]=1.0/9.0
t[6]=1.0/9.0
# 0 velocity weight is 4/9
t[0]=4.0/9.0
print('weights ',t)
#t[np.asarray([np.linalg.norm(ci)<1.1 for ci in c])] = 1./9.; t[0] = 4./9.
#print('weights ',t)
#
noslip = [lattice_vec.tolist().index((-lattice_vec[i]).tolist()) for i in range(q)] 
print('noslip onsite bounce back ',noslip)
i_set_neg_ux = np.arange(q)[np.asarray([ci[0]<0  for ci in lattice_vec])] # Unknown on right wall.
print('lattice vectors pointing to - x ',i_set_neg_ux)
i_set_0_ux = np.arange(q)[np.asarray([ci[0]==0 for ci in lattice_vec])] # Vertical middle.
print('lattice vectors with 0 x comp ',i_set_0_ux)
i_set_pos_ux = np.arange(q)[np.asarray([ci[0]>0  for ci in lattice_vec])] # Unknown on left wall.
print('lattice vectors pointing to + x i3',i_set_pos_ux)

#####
maxIter=20000#20000
print('run for time steps',maxIter)
# every sub
subIter=int(1.0/uBCLB) #int(maxIter/25)
#






# now call it
print('generating obstacles')

obstacle,binary=obstacle_def(r)


# Prep figure
fig, axs = plt.subplots(1,figsize=(8,5),dpi=120)
# for vector field plot, plot arrows every n_coarse sites
n_coarse=5
u=np.zeros((2,nx,ny))
ux_plot=np.transpose(u[0,::n_coarse,::n_coarse])*(1.0-np.transpose(binary[::n_coarse,::n_coarse]))
uy_plot=np.transpose(u[1,::n_coarse,::n_coarse])*(1.0-np.transpose(binary[::n_coarse,::n_coarse]))
#im_disc=axs.imshow(np.transpose(binary[::n_coarse,::n_coarse]), \
#                   cmap='Greens',alpha=np.transpose(binary[::n_coarse,::n_coarse]))
im_disc=axs.imshow(np.transpose(binary), \
                   cmap='Greens',alpha=np.transpose(binary))
#
print('shape of quiver plotting ',u[0,::n_coarse,::n_coarse].shape)
nx_plot,ny_plot=u[0,::n_coarse,::n_coarse].shape
x_plot=np.zeros((nx_plot,ny_plot))
y_plot=np.zeros((nx_plot,ny_plot))
for i in range(0,nx_plot):
    for j in range(0,ny_plot):
        x_plot[i,j]=float(i)*n_coarse
        y_plot[i,j]=float(j)*n_coarse
# NB making scale number bigger makes arrows smaller!
im=axs.quiver(np.transpose(x_plot),np.transpose(y_plot),ux_plot,uy_plot,scale=5)
#
vorticity=np.zeros((nx,ny))
im_vorticity=axs.imshow(np.transpose(vorticity),cmap='seismic', \
                      alpha=np.transpose(1.0-binary),vmin=-0.05,vmax=0.05)
#
axs.get_xaxis().set_visible(False)
axs.get_yaxis().set_visible(False)
axs.set_aspect('equal')
# set True to display to screen, False to save set of images to make movie
# NB saving images and displaying can be done but it may be that
# rescaling window alters images saved which if done during run can mess up
# movie made from theses images
display_not_save_images=True
#display_not_save_images=False
##
# set LHS BC
y=np.linspace(0,ny,ny)
sym_break_perturb=1.0e-1
ux_x0BC=uBCLB*(1.0+sym_break_perturb*np.exp(-(y-ny/3)**2))
if(abs(sym_break_perturb) > 1.0e-16): print('perturbing lh BC to break symmetry')


# initial velocities, 0 along y, mostly uBCLB along x
u=np.zeros((2,nx,ny))
u[0,:,:]=uBCLB
# lhs (ie inlet) set
u[0,0,:]=ux_x0BC
# initial conditions are equilibrium for rho = 1 and vel
rho_initial=1.0
feq = equilibrium(rho_initial,u)
# copy equilibrium f into initial f
fin = feq.copy()
###### Main time loop ##########################################################
#two_thirds=2.0/3.0
#one_sixth=1.0//6.0
for tstep in range(maxIter):
# Right wall: outflow condition. Set LB f in last two row on RHS to be equal
# i_set_neg_ux are the 3 vectors with negative x components
# in previous step PBCs will have set these values to values from LH, which need to be
# overwritten
    fin[i_set_neg_ux,-1,:] = fin[i_set_neg_ux,-2,:] 
# sum the 9 components of fin on each lattice site to get density rho at each lattice site
    rho = np.sum(fin,axis=0)         
# Calculate velocity.
#    u = np.dot(c.transpose(), fin.transpose((1,0,2)))/rho
# here for using numpy's einsum, convention is: i is over q lattice vectors,
# j=0 (x) and  1 (y) and m and n are over x and y lattice sites, respectively
    u=np.einsum('imn,ij',fin,lattice_vec,optimize=True)/rho
#    print(u.shape)
    '''
    Zou-He BC at left wall
    '''
# Left wall (ie. at x=0): 1st impose velocity
    u[0,0,:] = ux_x0BC
    u[1,0,:] = 0.0
# after streaming etc rho on left wall is wrong as np.roll applied PBCs there so fin
# with is from i_set_pos_ux are wrong as they rolled round from right hand side
# Zou and He (1997) showed that could compute rho starting from eqs for rho and u and
# eliminating f with i in i_set_pos_ux to get equation:
# compute density at left wall, ie for rho[0,:] overwrite values obtained by summing old fin above
    rho[0,:] = 1./(1.-u[0,0,:]) * \
          (np.sum(fin[i_set_0_ux,0,:],axis=0)+2.0*np.sum(fin[i_set_neg_ux,0,:],axis=0))
    '''
    now compute equilibrium fs, for all lattice including lh wal where we know have correct rho
    '''
    feq = equilibrium(rho,u)
# Left wall: Zou/He boundary condition
# Note that we don't know the fin for is with positive u_x for lh wall sites
# Zou-He said to apply bounce back but only to non-equil bit, eg apply
# bounce back to f_6 - feq_6 etc so here
# bounce back means f_6 - feq_6 = f_3 - feq_3, f_7 - feq_7 = f_4 - feq_4, f_8 - feq_8 = f_5 - feq_5
    fin[i_set_pos_ux,0,:] = fin[i_set_neg_ux,0,:] + feq[i_set_pos_ux,0,:] - feq[i_set_neg_ux,0,:]
# NB can also write fin explictly in terms of know fins, rho and u, see 1997 paper of Zou-He
# that is equivalent to just applying boundeback to non-equil parts of the three f_i, i=6,7,8 
# as we do here
# NB2 In original Palabos code last feq in line above is fin, so fin is just set equal to
# feq at given rho and u, i.e., Palabos code does not quite implement Zou-He BC
# but this difference seems to have very little affect, get same results with it as with
# normal Zou-He BC
# Collision step
    fout = fin - omega * (fin - feq)
# BCs at obstacles: flip all vectors for all cells in obstacle
# eg f for [0 -1] becomes f for [0 1]
    for i in range(q): 
        fout[i,obstacle] = fin[noslip[i],obstacle]
# Streaming step with PBCs implemented using numpy#'s roll command
    for i in range(q):
        fin[i,:,:] = np.roll(np.roll(fout[i,:,:],lattice_vec[i,0],axis=0),lattice_vec[i,1],axis=1)
# write out stuff
    if(tstep > 1 and tstep%subIter == 0):
        animation_plot(display_not_save_images)
        mean_ux=np.mean(u[0,:,:])#/uBCLB)
        min_ux=np.amin(u[0,:,:])#/uBCLB)
        max_ux=np.amax(u[0,:,:])#/uBCLB)
        print('time step ',tstep,' time/(2*r/uBCLB)',round(tstep/T_d,2), \
              'mean, min and max u_x ',round(mean_ux,6),round(min_ux,6),round(max_ux,6))
#
tend = time.time()
print('runtime = ',round((tend-tstart)/60.0,2),' minutes')
