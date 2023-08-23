from simulation import Simulation

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import constants
from scipy.signal import savgol_filter, get_window

class Spectra(Simulation):
    c = constants.value('speed of light in vacuum')
    def __init__(self, 
                 vacf_file = '', 
                 dipole_file = '',
                 dt = 0.25):

        super().__init__(vacf_file = '', dipole_file = '')
        self.vacf_file = vacf_file
        self.dipole_file = dipole_file
        self.dt = dt


    def HA(self, V, T):
        """
        Apply rescaling of the relative line intensities according to the
        harmonic approximation. For further inside see:
        https://doi.org/10.1063/1.1774986

        Parameters:
        V (list):       array of wavenumbers in cm^-1
        T (float):      Temperature in K

        Returns:
        Q (list):       SCaling function
        """
        conv = 0.123985e-3  # scaling factor: cm^-1 in eV (1/8065.5)
        kB = 8.617333e-5    # Boltzmann constant in eV/K
        hbar = 6.582119569e-16  # Planck's constant in eV*s
        beta = 1 / (kB * T)
        Q = (beta * hbar * np.array(V) * conv) / (1 - np.exp(-1*beta * hbar * V * conv))

        return Q


    def get_suffix(self, s):
        try:
            suffix = s.split('/')[-1].split('.')[-2]
        except IndexError:
            suffix = s.split('/')[-1].split('.')[-1]
        
        return suffix


    def get_mode_spectrum(self, t, acf):
        """
        Calculate the fourier transform of some autocorrelated data 'acf',
        i.e. the mode spectrum with respect to d. Also applies Blackmann-
        Harris window function.

        Parameters:
        t (list):       time in fs
        acf (list):     autocorrelated data

        Returns:
        v (list):       wavenumbers in cm^-1
        fft_acf (list): Fourier transform of acf
        """
        N = t.size
        # normalize vacf
        norm_acf = acf / max(acf)
        # fourier transform and apply window function
        fft_acf = np.abs(np.fft.rfft(norm_acf * get_window(
                    'blackmanharris', N)))
        # wavenumber in cm^-1
        v = 0.01 * np.fft.rfftfreq(N, d = self.dt) * 1e+15 / self.c

        return v, fft_acf


    def velocity(self, plot = False):
        """
        Computes the mode spectrum from the velocity autocorrelation function.
        Units in the 'vacf_file' are expected to be in LAMMPS real units
        """
        try:
            if not self.vacf_file:
                raise FileNotFoundError('No vacf_file provided for VACF spectrum calculation...')
        except FileNotFoundError as e:
            exit(e)

        ### load data
        data = np.genfromtxt(self.vacf_file,
                            skip_header = 1)
        t = data[:, 0]              # time in fs
        vacf = data[:, -1]          # vacf in (Angstr/fs)^2
        
        ### get the file name's suffix
        suffix = self.get_suffix(self.vacf_file)

        ### calculate the mode spectrum
        v, fft_vacf = self.get_mode_spectrum(t, vacf)

        ### filtering
        # apply savgol filter: 
        # 'window_size' and 'polyorder' are hardcored here... maybe find
        # some way to automatically determine them 
        fft_vacf_filtered = savgol_filter(fft_vacf, 500+1, 5, mode='mirror')
        fft_vacf = fft_vacf_filtered

        ### save to file
        df_out = pd.DataFrame({'v [cm^-1]': v, 'fft_vacf': fft_vacf})
        df_out.to_csv(f'{suffix}.dat', sep = '\t',
                     float_format = '%.4f', index = False)
        
        if plot:
            ### plot the data
            fig, (ax1, ax2) = plt.subplots(2, 1)
            # upper plot: VACF
            ax1.plot(t, vacf, '-', lw = 0.5)
            ax1.set_xlabel(r'Time $t$ [fs]')
            ax1.set_ylabel(r'VACF [$\AA^2$/fs$^2$]')
            # lower plot: spectrum
            ax2.plot(v, fft_vacf, '-', lw = 1.0)
            ax2.set_xlim((0, 10000))    # >10000 cm^-1 usually not interesting...
            ax2.set_xlabel(r'Wavenumber $\tilde{\nu}$ [cm$^{-1}$]')
            ax2.set_ylabel(r'|$\mathcal{F}$ (VACF)| [arb. units]')

            fig.suptitle(r'VACF & $\mathcal{F}$(VACF)')
            plt.tight_layout()
            plt.savefig(f'{suffix}.png', bbox_inches='tight', dpi=600)


    def dipole(self, plot = False, T = 300):
        """
        Computes the mode spectrum from the change in dipole moment. This
        is then effectively the IR spectrum. For intensity rescaling the
        harmonic approximation is used (https://doi.org/10.1063/1.1774986).
        Units in the 'dipole_file' are expected to be in LAMMPS real units
        """
        try:
            if not self.dipole_file:
                raise FileNotFoundError('No dipole_file provided for DACF spectrum calculation...')
        except FileNotFoundError as e:
            exit(e)

        ### load data
        data = np.genfromtxt(self.dipole_file,
                            skip_header = 1)
        t = data[:-1, 0]               # time in fs
        d_arr = data[:, -1]          # dipole moment in e*Angstroms
        
        ### get the file name's suffix
        suffix = self.get_suffix(self.dipole_file)

        ### autocorrelate the dipole moment's derivative
        d_diff = np.diff(d_arr)
        dacf = np.correlate(d_diff, d_diff, mode='full')[len(d_diff)-1:]

        ### calculate the mode spectrum
        v, fft_dacf = self.get_mode_spectrum(t, dacf)
        
        ### rescale the lineshape (see HA method for details)
        Q = np.array(self.HA(v, T))
        ir = np.multiply(fft_dacf, Q)

        ### filtering
        # apply savgol filter: 
        # 'window_size' and 'polyorder' are hardcored here... maybe find
        # some way to automatically determine them 
        ir_filtered = savgol_filter(fft_dacf, 500+1, 5, mode='mirror')
        ir = ir_filtered

        ### save to file
        df_out = pd.DataFrame({'v [cm^-1]': v, 'fft_dacf': fft_dacf})
        df_out.to_csv(f'{suffix}.dat', sep = '\t',
                     float_format = '%.4f', index = False)

        if plot:
            ### plot the data
            fig, (ax1, ax2) = plt.subplots(2, 1)
            # upper plot: DACF
            ax1.plot(t, dacf, '-', lw = 0.5)
            ax1.set_xlabel(r'Time $t$ [fs]')
            ax1.set_ylabel(r'ACF($\dot{d}$) [$e^2\AA^2/fs^2$]')
            # lower plot: spectrum
            ax2.plot(v, fft_dacf, '-', lw = 1.0)
            ax2.set_xlim((0, 10000))    # >10000 cm^-1 usually not interesting...
            ax2.set_xlabel(r'Wavenumber $\tilde{\nu}$ [cm$^{-1}$]')
            ax2.set_ylabel(r'|$\mathcal{F}$ (ACF)| [arb. units]')

            fig.suptitle(r'DACF & $\mathcal{F}$(DACF)')
            plt.tight_layout()
            plt.savefig(f'{suffix}.png', bbox_inches='tight', dpi=600)


