import numpy as np

def particle_vel():
    #function that returns particle velcoity for a given time
    #uses drag+gravity in Maxey & Riley equation
    #implements forward difference finite difference method
    up = np.zeros(len(time))
    for it in range(len(time)):
        up[it] = (del_t/m_p) * (-6*np.pi*mu*0.5*dp*(uf-up[it-1]) +((m_p-m_f)*g))

    return up

def particle_history():
    #function that returns particle velcoity for a given time
    #uses drag+gravity in Maxey & Riley equation
    #implements forward difference finite difference method
    up = np.zeros(len(time))
    for it in range(len(time)):

        up[it] = (del_t/m_p) * (-6*np.pi*mu*0.5*dp*(uf-up[it-1]) +((m_p-m_f)*g) + F_b)

    return up



if __name__ =="__main__":
    #print("Booyeah")
    #diameter, reynold numb of particle
    dp = 1.0
    Re_p = 0.1
    #densities
    rho_p = 1.8
    rho_f = 1.0
    #gravity
    g = 1.0
    #fluid viscosities
    mu = np.sqrt(rho_f*rho_p*(1-rho_f/rho_p)*g*(dp**3) /(18*Re_p) )
    nu = mu/ rho_f
    #particle's volume
    Vp=(np.pi * dp**3) / 6.0
    #mass of particle
    m_p = rho_p*Vp 
    #displaced mass of fluid by the presence of particle.
    m_f =  rho_f*Vp
    #particle relaxation time.
    tau_p=rho_p*dp**2/(18*mu) 
    N_tot=256
    time=np.linspace(0,2*np.pi,N_tot)
    del_t = time[1]-time[0]
    #quiescent fluid
    uf =0
    #particle velocity calculated
    vel = particle_vel()
    print(range(len(time)))
    print(vel[-5:-1])

    
    
    



    
     
