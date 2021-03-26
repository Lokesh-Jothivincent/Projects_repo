import numpy as np

def particle_vel():
    #function that returns particle velcoity for a given time
    #uses drag+gravity in Maxey & Riley equation
    #implements forward difference finite difference method
    up = np.zeros(len(time))
    for it in range(1,len(time)):
        up[it] = (del_t/m_p) * (-6*np.pi*mu*0.5*dp*(uf[it-1]-up[it-1]) +((m_p-m_f)*g))

    return up

def particle_vel_history():
    #function that returns particle velcoity for a given time
    #uses drag+gravity in Maxey & Riley equation
    #impleemnts forward difference finite difference method
    up = np.zeros(len(time))
    g_t = np.zeros(len(time))
    C_b = 6 * (0.5*dp**2) * rho_f *np.sqrt(np.pi * nu)

    for it in range(1,len(time)):
        N=it-1
        start_ind = N- (N_win -1)
        if it == 1:
            F_h = (4/3)*C_b * g_t[N]*del_t

        elif N_win <= it:
            #for loop for calculating third term
            temp_sum = 0
            rev_ind = N-1
            for n in range(1,N-1):
                temp_sum += g_t[rev_ind] *(  ((n+(4/3)) / ((n+1)**1.5 +(n+1.5)*(n**0.5)) )\
                    + ((n-(4/3)) / ((n-1)**1.5 +(n-1.5)*(n**0.5) ))   )
                rev_ind -=1
            #history force term
            F_h = (4/3)*C_b * g_t[N]*del_t + C_b*g_t[0] * \
                ( ( N-(4/3) )/( (N-1)**1.5 +(N-1.5)* (N**0.5) )) \
                    +C_b * del_t * temp_sum

        else:
            #starting index
            start_ind = N- (N_win -1)
            #for loop for calculating third term
            temp_sum = 0
            #rev_index for calculating acceleration g_t
            rev_ind = N-1
            #third term summation
            for n in range(1,N-1):
                temp_sum += g_t[rev_ind] *(  ((n+(4/3)) / ((n+1)**1.5 +(n+1.5)*(n**0.5)) )\
                    + ((n-(4/3)) / ((n-1)**1.5 +(n-1.5)*(n**0.5) ))   )
                rev_ind -=1
            #history force
            F_h = (4/3)*C_b * g_t[N]*del_t + C_b*g_t[start_ind] * \
                ( ( N-(4/3) )/( (N-1)**1.5 +(N-1.5)* (N**0.5) )) \
                    +C_b * del_t * temp_sum

        #particle velocity
        up[it] = (del_t/m_p) * (-6*np.pi*mu*0.5*dp*(uf[it-1]-up[it-1]) +((m_p-m_f)*g ) +F_h)
        #acceleration
        g_t[it] = ((uf[it] - uf[it-1])/del_t) - ((up[it] - up[it-1])/del_t)

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
    uf =np.zeros(len(time))
    #particle velocity calculated
    vel = particle_vel()
    print(vel[-5:-1])
    N_win = 5 #python index from 0 to N-1 5points here from 0 to N-1
    vel_history = particle_vel_history()
    print(vel_history[-5:-1])
    

    
    
    



    
     
