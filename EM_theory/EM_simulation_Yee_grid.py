import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
import pdb 

sstyle = "seaborn-poster"
plt.style.use(sstyle)
plt.rc('font',family = 'serif')

ε0 = spc.epsilon_0
μ0 = spc.mu_0
c  = spc.speed_of_light
qe = spc.elementary_charge
me = spc.electron_mass


def EM_wave_1D(Nt, courant = 1, counter_prop=False):
    """
    This function aims to simulate a 1D electromagnetic wave. 
    The function can also simulate counter propagating waves 
    such that there is no initial magnetic fiel.  

    Nt: number of steps taken as time. DType: Int.
    courant: courant number, depending on the factor, is the
             stability and efficiency of algorithm.
             Defaul: 1, usually is good when <= 1, for stability,
    counter_prop: Defautl: False.
    """
    σ  = courant

    Δz = 1.0; Nz = 500; 
    Δt = σ*Δz/c

    z_fp = np.arange(Nz)*Δz # spatial axis for full-points
    z_mp = z_fp[:-1] + 0.5*Δz  # spatial axis for mid-points

    kz = 2*np.pi/(16*Δz)    # wavenumber of Gaussian wavepacket 
    dk = 10/kz              # width of Gaussian wavepacket
    ff = c*kz/(2*np.pi)     # freq of Gaussian wavepacket
    # Initial conditions and setting up the E-field
    E0  = np.cos(kz * z_fp - 0.0 * Δt * ff * 2 * np.pi)*\
            np.exp( - ( (z_fp - 250 * Δz) / dk ) ** 2)
    E0[0], E0[-1] = 0, 0    # reflecting boundaries (?)          
    E  = E0.copy()  
    Et = np.zeros((Nt,Nz)) # all-times for E-field
    # Initial conditions and setting up the H-field       
    if counter_prop == False:
        H0 = np.sqrt(ε0/μ0)*np.cos( kz*z_mp - 0.5*Δt*ff*2*np.pi)*\
                np.exp( - ( (z_mp - 250.5 * Δz) / dk )**2) 
        H = H0.copy()
        numsteps = 10
    else:
        numsteps = 1
        H = np.zeros(Nz-1)

    def vacuumstep1d( E, H, numsteps):
        for i in range(numsteps):
            E[1:-1] += (Δt/(ε0*Δz))*( H[:-1] - H[1:] )
            H += (Δt/(μ0*Δz))*( E[:-1] - E[1:] )
        return 

    for n in range(Nt):
        E[1:-1] += (Δt/(ε0*Δz))*( H[:-1] - H[1:] )
        H += (Δt/(μ0*Δz))*( E[:-1] - E[1:] )
        Et[n,:] = E
        if (n%10) == 9:
            plt.clf()        
            plt.subplot(1,2,1)
            plt.plot(z_fp,E0, c="k", label="Initial Electric Field")
            plt.plot(z_fp,E, c="r", label="Electric Field")
            plt.title(f"t = {n}")
            plt.xlabel("Distance $z$")
            plt.ylabel("Amplitude $E_x$")
            plt.ylim(-1.5, 1.5)
            plt.legend(loc=1)
            plt.subplot(1,2,2)
            plt.plot(z_mp,H, c="b", label="Magnetic Field")
            plt.title(f"t = {n}")
            plt.xlabel("Distance $z$")
            plt.ylabel("Amplitude $H_y$")
            plt.ylim(-0.0027, 0.0027)
            plt.legend(loc=1)
            plt.tight_layout()
            plt.pause( 0.00001) 

    # plot final
    plt.clf()
    plt.plot( E, label="Final E-field")
    plt.plot( -E0[::-1], c="k", label="Initial E-field")   # reversed in time and space
    plt.xlabel("Distance $z$")
    plt.ylabel("Amplitude")
    plt.legend()
    plt.show()

    """
    # plot dispersion relation
    plt.clf()
    plt.imshow(abs(np.fft.fft2(Et)))
    plt.title("Dispersion Relation")
    plt.colorbar()
    plt.show() 
    """

EM_wave_1D(Nt = 1000, courant=0.9)

def third_and_fourth_problem():
    dt_ = float(input("set dt-factor for magic step: "))
    ω_pe_ = float(input("set ω_pe-factor: "))

    Δz = 1.0; Nz = 1000; S = 1
    Δt_magic = S*Δz/c
    Δt = dt_*Δt_magic
    Nt = round(Nz/(c*Δt))

    print(f"Number of total step-time: {Nt}")
     
    z  = np.arange(Nz)*Δz # spatial axis

    """
    n0 = 10e18/(0.01**3) # initial electron density: m^-3
    r  = lambda i: np.sqrt(3)*i*Δz 
    n  = lambda i: n0*(1 - (r(i)/0.5)**2) 
    ω_pe = lambda i: np.sqrt( n(i)*qe*qe/(me*ε0) ) #  plasma frequency
    """

    kz = 2*np.pi/(16*Δz); dk = 10/kz # wavenumber and width of Gaussian wave-packet
    ff = c*kz/(2*np.pi) # freq of Gaussian wave-packet
    ω_pe = np.array([0]*(Nz//2) +  [ω_pe_*2*np.pi*ff]*(Nz//2)) #  plasma frequency

    E = np.cos(kz * z - 0.0 * Δt * ff * 2 * np.pi)*\
         np.exp(( - ( (z - 250 * Δz) / dk ) ** 2)) # Initial condition on E
    H = np.zeros(Nz) # Also, is it better to have it as E/c instead?
    Jx = np.zeros(Nz) 

    def vacuumstep1d( E, H, J, i):
        J_ = J[i] + Δt*ε0*ω_pe[i]*ω_pe[i]*E[i]
        E_ = E[i+1] - (Δt/(ε0*Δz))*(H[i+1] - H[i]) - Δt*J[i]
        H_ = H[i] - (Δt/(μ0*Δz))*(E[i+1] - E[i])
        return E_,H_,J_

    for n in range(Nt):
        for i in range(Nz-1):
            E[i+1], H[i], Jx[i] = vacuumstep1d(E,H,Jx,i)
        E[0] = 0; E[-1] = 0; Jx[0] = 0; Jx[-1] = 0     
        plt.clf()        
        plt.subplot(1,2,1)
        plt.plot(z,E, c="r", label="Electric Field")
        plt.title(f"t = {n}")
        plt.xlabel("Distance $x$")
        plt.ylabel("Amplitude $E_x$")
        plt.ylim(-1, 1)
        plt.legend(loc=1)
        plt.subplot(1,2,2)
        plt.plot(z[:-1] + 0.5, H[:-1], c="b", label="Magnetic Field")
        plt.title(f"t = {n}")
        plt.xlabel("Distance $x$")
        plt.ylabel("Amplitude $H_y$")
        plt.ylim(-0.0027, 0.0027)
        plt.legend(loc=1)
        plt.tight_layout()
        plt.pause( 0.0001) 

"""p = input("which problem do you want to run? [1/2/3]: ")
if p == "1":
    cp = bool(input("Case of counter propagation? True/False:"))
    first_problem(counter_prop=cp)
elif p =="2":
    second_problem()
elif p =="3":
    third_and_fourth_problem()    
"""