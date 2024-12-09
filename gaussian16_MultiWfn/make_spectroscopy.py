import numpy as np
from matplotlib import pyplot as plt
import subprocess as sp

def get_globals():
    global OUT_FILE, PLOTS_DIR
    OUT_FILE = "geometry.out" # Gaussian16 output file

    PLOTS_DIR = "PLOTS_DATA"
    sp.call(f"mkdir -p {PLOTS_DIR}", shell=True)

def read_data():
    global NSTATES, EL_DIP, MAG_DIP, ENERGY, OSC_STR, ROT_STR, EM_ANG
    NSTATES = int( sp.check_output("grep 'Excited State' %s | tail -n 1 | awk '{print $3}'" % (OUT_FILE), shell=True).decode().split(":")[0])
    EL_DIP  = np.zeros((NSTATES, 3))
    MAG_DIP = np.zeros((NSTATES, 3))
    ENERGY  = np.zeros((NSTATES))
    OSC_STR = np.zeros((NSTATES))
    ROT_STR = np.zeros((NSTATES))
    EM_ANG  = np.zeros((NSTATES))

    sp.call( "grep 'transition electric dipole' %s -A %d | tail -n %d | awk '{print $2, $3, $4}' > temp.dat" % (OUT_FILE, NSTATES+1, NSTATES), shell=True )
    EL_DIP += np.loadtxt("temp.dat") # in a.u.
    sp.call( "grep 'transition magnetic dipole' %s -A %d | tail -n %d | awk '{print $2, $3, $4}' > temp.dat" % (OUT_FILE, NSTATES+1, NSTATES), shell=True )
    MAG_DIP += np.loadtxt("temp.dat") # in a.u.
    for state in range(NSTATES):
        ENERGY[state]  = float( sp.check_output( "grep 'Excited State' %s | tail -n %d | head -n 1 | awk '{print $5}'" % (OUT_FILE,NSTATES-state), shell=True ).decode() )
        OSC_STR[state] = float( sp.check_output( "grep 'Excited State' %s | tail -n %d | head -n 1 | awk '{print $9}'" % (OUT_FILE,NSTATES-state), shell=True ).decode().split("=")[1] )
    #ENERGY /= 27.2114 # eV --> a.u.
    sp.call( "grep 'R(velocity)' %s -A %d | tail -n %d | awk '{print $5}' > temp.dat" % (OUT_FILE, NSTATES, NSTATES), shell=True )
    ROT_STR += np.loadtxt("temp.dat") # in cgs (10**-40 erg-esu-cm/Gauss)
    sp.call( "grep 'R(velocity)' %s -A %d | tail -n %d | awk '{print $6}' > temp.dat" % (OUT_FILE, NSTATES, NSTATES), shell=True )
    EM_ANG += np.loadtxt("temp.dat") # in cgs (10**-40 erg-esu-cm/Gauss)
    sp.call("rm temp.dat", shell=True)

    np.savetxt(f"{PLOTS_DIR}/ENERGY_OSC_ROT.dat", np.c_[ENERGY, OSC_STR, ROT_STR], fmt="%12.6f", header="Energy(eV) OSC_STR ROT_STR")
    np.savetxt(f"{PLOTS_DIR}/EL_DIP_MAG_DIP_ANGLE.dat", np.c_[EL_DIP, MAG_DIP, EM_ANG], fmt="%12.6f", header="EL_DIP_x EL_DIP_y EL_DIP_z MAG_DIP_x MAG_DIP_y MAG_DIP_z EM_ANGLE")



def plot_data():

    NPTS  = 5000
    EMIN  = np.min(ENERGY) - 0.1 # eV
    EMAX  = np.max(ENERGY) + 0.1 # eV
    EGRID = np.linspace(EMIN, EMAX, NPTS)
    SIGMA = 0.05 # eV
    ABS_G   = np.zeros((NPTS))
    CD_G    = np.zeros((NPTS))
    ABS_L   = np.zeros((NPTS))
    CD_L    = np.zeros((NPTS))
    for pt in range(NPTS):
        ABS_G[pt] += np.sum( OSC_STR[:] * np.exp( -1.0 * (EGRID[pt] - ENERGY[:])**2 / (2.0 * SIGMA**2) ) )
        CD_G[pt]  += np.sum( ROT_STR[:] * np.exp( -1.0 * (EGRID[pt] - ENERGY[:])**2 / (2.0 * SIGMA**2) ) )
        ABS_L[pt] += np.sum( OSC_STR[:] *  SIGMA**2/4/( (EGRID[pt] - ENERGY[:])**2  + (SIGMA/2)**2) )
        CD_L[pt]  += np.sum( ROT_STR[:] *  SIGMA**2/4/( (EGRID[pt] - ENERGY[:])**2  + (SIGMA/2)**2) )

    np.savetxt(f"{PLOTS_DIR}/ABS_CD_spectra.dat", np.c_[EGRID, ABS_G, CD_G, ABS_L, CD_L], fmt="%12.6f", header="Energy(eV) Abs_G CD_G Abs_L CD_L")
    np.savetxt(f"{PLOTS_DIR}/EM_ANGLE.dat", np.c_[ENERGY, EM_ANG], fmt="%12.6f", header="Energy(eV) EM_Angle(deg)")


    # Absorption on left vetical axis and rotary strength on right vertical axis
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(EGRID, ABS_G, c="black", label="Absorption")
    #ax1.stem(ENERGY, OSC_STR, linefmt="black", markerfmt="black", basefmt=" ")
    ax1.scatter(ENERGY, OSC_STR, edgecolors='black', facecolors='none', marker="o", linewidth=2, s=40)
    ax2.plot(EGRID, CD_G, c="red", label="CD")
    #ax2.stem(ENERGY, ROT_STR, linefmt="red", markerfmt="red", basefmt=" ")
    ax2.scatter(ENERGY, ROT_STR, edgecolors='red', facecolors='none', marker="o", linewidth=2, s=40)
    #fig.legend(bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)
    ax1.set_xlabel("Energy (eV)", fontsize=15)
    ax1.set_ylabel("Absorption (Osc. Str.)", fontsize=15)
    ax2.set_ylabel("Circular Dichroism (Rot. Str.)", fontsize=15)
    plt.xlim(EMIN, EMAX)
    ax1.set_ylim(-np.max(ABS_G)*1.1,np.max(ABS_G)*1.1)
    plt.tight_layout()
    plt.savefig(f"{PLOTS_DIR}/ABS_CD_spectra_GAUSSIAN.jpg", dpi=300)
    plt.clf()
    plt.close()

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(EGRID, ABS_L, c="black", label="Absorption")
    #ax1.stem(ENERGY, OSC_STR, linefmt="black", markerfmt="black", basefmt=" ")
    ax1.scatter(ENERGY, OSC_STR, edgecolors='black', facecolors='none', marker="o", linewidth=2, s=40)
    ax2.plot(EGRID, CD_L, c="red", label="CD")
    #ax2.stem(ENERGY, ROT_STR, linefmt="red", markerfmt="red", basefmt=" ")
    ax2.scatter(ENERGY, ROT_STR, edgecolors='red', facecolors='none', marker="o", linewidth=2, s=40)
    #fig.legend(bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)
    ax1.set_xlabel("Energy (eV)", fontsize=15)
    ax1.set_ylabel("Absorption (Osc. Str.)", fontsize=15)
    ax2.set_ylabel("Circular Dichroism (Rot. Str.)", fontsize=15)
    plt.xlim(EMIN, EMAX)
    ax1.set_ylim(-np.max(ABS_G)*1.1,np.max(ABS_G)*1.1)
    plt.tight_layout()
    plt.savefig(f"{PLOTS_DIR}/ABS_CD_spectra_LORNTZIAN.jpg", dpi=300)
    plt.clf()
    plt.close()

    # Plot the angles
    fig, ax1 = plt.subplots()
    ax1.plot(EGRID, ABS_L/np.max(ABS_L)*np.max(EM_ANG), c="black", label="Abs")
    ax1.plot(EGRID, CD_L/np.max(np.abs(CD_L))*np.max(EM_ANG), c="red", label="CD")
    ax1.plot(ENERGY, EM_ANG, "-o", c="green", lw=3, ms=10, label="EM Angle")
    ax1.set_xlabel("Energy (eV)", fontsize=15)
    ax1.set_ylabel("EM Angle (deg)", fontsize=15)
    fig.legend(bbox_to_anchor=(1,1), bbox_transform=ax1.transAxes)
    plt.tight_layout()
    plt.savefig(f"{PLOTS_DIR}/EM_ANGLE.jpg", dpi=300)
    plt.clf()
    plt.close()




def main():
    get_globals()
    read_data()
    plot_data()

if ( __name__ == "__main__"):
    main()