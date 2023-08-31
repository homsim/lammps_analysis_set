from simulation import Simulation

from ovito.io import import_file
from ovito.data import BondsEnumerator
from ovito.modifiers import CreateBondsModifier, CoordinationAnalysisModifier
from ovito.modifiers import ExpressionSelectionModifier
import numpy as np
from numpy.linalg import norm
import pandas as pd
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

class Temperatures(Simulation):
    """
    This class holds timeseries for center-of-mass, vibrational and rotational
    velocity components, as well as their respective temperature.
    Maybe also the distributions (dont know how to handle time here yet)...
    """
    kB = 1.380649e-23           # Boltzmann's constant in J/K
    N_A = 6.02214076e+23        # Avogadro's constant
    conv = 1.6605390671738466e-17               # conversion:
    #conv = (1e-3/N_A)*((1e-10)**2/(1e-15)**2)  # (g*Angstrom^2)/(mole*fs^2) to J 
                                                
    def __init__(self, 
                 traj_file = '',
                 bond_file = '',
                 dt = 0.25,
                 types = []):

        super().__init__(traj_file = '', bond_file = '')
        self.traj_file = traj_file
        self.bond_file = bond_file
        self.dt = dt
        self.types = types

    def get_modifiers(self, data):
        """
        This method defines the ovito modifiers that determine which particle
        is identified as a gas molecule. It also creates the bonds needed
        to assign the bond properties to.
        """
        modifiers = []

        # detect if a ReaxFF bond file is given and if yes, take the bonds 
        # from it. If not, create bonds according atoms' vdW radii
        if not self.bond_file:
            ### get vdw_radius to use in the bond cutoff
            r_vdw_1 = data.particles.particle_types.types[self.types[0] - 1].vdw_radius
            r_vdw_2 = data.particles.particle_types.types[self.types[1] - 1].vdw_radius
            # arbitrarily define the bond cutoff as the mean of bond vdw radii
            # time 1.2 -> 1.2 / 2 = 0.6 (consistent with Lorentz-Berthelot rules)
            bond_cutoff = 0.6 * (r_vdw_1 + r_vdw_2) 

            # create_bonds
            create_bonds_mod = CreateBondsModifier(
                mode=CreateBondsModifier.Mode.Pairwise)
            create_bonds_mod.set_pairwise_cutoff(
                self.types[0], self.types[1], bond_cutoff) 
            modifiers.append(create_bonds_mod)
        else:
            load_traj_mod = LoadTrajectoryModifier()
            load_traj_mod.source.load(self.bond_file)
            modifiers.append(load_traj_mod)

        # identify diatomic molecules in the gas phase via coordination
        # and select them
        coord_an_mod = CoordinationAnalysisModifier()
        modifiers.append(coord_an_mod)
        expr_sel_mod = ExpressionSelectionModifier(
            expression='Coordination<2')
        modifiers.append(expr_sel_mod)

        return modifiers

    def create_particle_props(self, data):
        """
        Create the particle properties needed for the mode velocities.
        """
        N_part = data.particles.count
        null_data = np.zeros((N_part, 3))     # to place in the new props temporarily
        data.particles_.create_property('v_com', data = null_data)
        data.particles_.create_property('v_vib', data = null_data)
        data.particles_.create_property('v_rot', data = null_data)
        data.particles_.create_property('R', data = null_data)
        data.particles_.create_property('d', data = null_data)

        return 0

    def create_bond_props(self, data):
        """
        Create the bond properties needed for the mode temperatures.
        """
        N_bond = data.particles.bonds.count
        null_data = np.zeros(N_bond)
        data.particles_.bonds_.create_property('T_com', data = null_data)
        data.particles_.bonds_.create_property('T_vib', data = null_data)
        data.particles_.bonds_.create_property('T_rot', data = null_data)

        return 0

    def iterate_frames(self, plot = False):
        """
        Iterate over all frames in the traj_files and calculate mode
        velocities and temperatures.
        """
        try:
            if not self.traj_file:
                raise FileNotFoundError('No traj_file provided for temperature calculation...')
        except FileNotFoundError as e:
            exit(e)

        pipeline = import_file(self.traj_file)

        data = pipeline.compute()

        mods = self.get_modifiers(data)
        for mod in mods:
            pipeline.modifiers.append(mod)

        t_arr = []
        T_com_arr = []
        T_vib_arr = []
        T_rot_arr = []
        for frame in range(pipeline.source.num_frames):
            data = pipeline.compute(frame)
            
            t_arr.append(data.attributes['Timestep']*self.dt/1000)

            ### create particle properties for the velocities
            self.create_particle_props(data)

            self.get_mode_velocities(data)
            self.get_mode_temperatures(data)
            
            # filter zeros, i.e. from non-gas phase particles/bonds
            #self.v_com = [v for v in data.particles['v_com'] if any(v)]
            #self.v_vib = [v for v in data.particles['v_vib'] if any(v)]
            #self.v_rot = [v for v in data.particles['v_rot'] if any(v)]
            T_com = np.mean([T for T in data.particles.bonds['T_com'] if T])
            T_vib = np.mean([T for T in data.particles.bonds['T_vib'] if T])
            T_rot = np.mean([T for T in data.particles.bonds['T_rot'] if T])

            T_com_arr.append(T_com)
            T_vib_arr.append(T_vib)
            T_rot_arr.append(T_rot)

        df_out = pd.DataFrame({'t [ps]': t_arr,
                               'T_com': T_com_arr,
                               'T_vib': T_vib_arr,
                               'T_rot': T_rot_arr})
        df_out.to_csv('Mode_temp_gas.dat', sep = '\t',
                      float_format = '%.2f', index = False)

        if plot:
            fig, ax = plt.subplots(1, 1, figsize = (6.4, 4.8))
            ax.plot(t_arr, T_com_arr, label = '$T_{COM}$')
            ax.plot(t_arr, T_vib_arr, label = '$T_{vib}$')
            ax.plot(t_arr, T_rot_arr, label = '$T_{rot}$')

            ax.set_xlabel('Time $t$ [ps]')
            ax.set_ylabel('Mode temp. $T$ [K]')

            ax.legend(loc = 'best')
            plt.tight_layout()
            plt.savefig('Mode_temp_gas.png', dpi = 600, bbox_inches = 'tight')

        return 0

    def get_mode_temperatures(self, data):
        ### get relevant properties
        v_coms = data.particles['v_com']
        v_vib = data.particles['v_vib']
        v_rot = data.particles['v_rot']
        R = data.particles['R']
        d = data.particles['d']
        ptypes = data.particles['Particle Type']

        ### create bond property for the temperatures
        self.create_bond_props(data)

        sel = data.particles['Selection']

        ### iterate over the created bonds' topology
        bond_topology = data.particles.bonds.topology
        bonds_enum = BondsEnumerator(data.particles.bonds)
        for a,b in data.particles.bonds.topology:
            if sel[a] == 1 and sel[b] == 1:
                # analyse only the type combinations specified in 'types'
                type_comb = sorted([ptypes[a], ptypes[b]])
                if type_comb == sorted(self.types):
                    # world's most complicated way to get particle masses...
                    m_a = data.particles.particle_types.types[ptypes[a] - 1].mass
                    m_b = data.particles.particle_types.types[ptypes[b] - 1].mass

                    v_com = v_coms[a]
                    v_vib_a = v_vib[a]
                    v_vib_b = v_vib[b]
                    v_rot_a = v_rot[a]
                    v_rot_b = v_rot[b]
                    
                    d_ab = d[a]
                    d_ba = d[b]
                    R_a = R[a]
                    R_b = R[b]

                    # get the index of the bond that connects a and b only
                    bonds_of_a = bonds_enum.bonds_of_particle(a)
                    bonds_of_b = bonds_enum.bonds_of_particle(b)
                    bond_index = next(iter(set(bonds_of_a).intersection(bonds_of_b)))
                    
                    # calculate the temperatures
                    E_com = .5 * (m_a + m_b) * norm(v_com)**2
                    T_com = 0.66666666 * E_com * self.conv / self.kB
                    data.particles_.bonds['T_com'][bond_index] = T_com 

                    E_vib = .5 * (m_a * norm(v_vib_a)**2 
                                 + m_b * norm(v_vib_b)**2)
                    T_vib = 2 * E_vib * self.conv / self.kB
                    data.particles_.bonds['T_vib'][bond_index] = T_vib

                    I = m_a * norm(R_a)**2 + m_b * norm(R_b)**2
                    omega_a = 1/norm(R_a)**2 * np.cross(R_a, v_rot_a)
                    E_rot = .5 * I * norm(omega_a)**2
                    T_rot = E_rot * self.conv / self.kB
                    data.particles_.bonds['T_rot'][bond_index] = T_rot
        
        return 0

    def get_mode_velocities(self, data):
        ### get relevant properties
        positions = data.particles['Position']
        velocities = data.particles['Velocity']
        ptypes = data.particles['Particle Type']

        sel = data.particles['Selection']

        ### iterate over the created bonds' topology
        bond_topology = data.particles.bonds.topology
        bonds_enum = BondsEnumerator(data.particles.bonds)
        for a,b in data.particles.bonds.topology:
            if sel[a] == 1 and sel[b] == 1:
                # analyse only the type combinations specified in 'types'
                type_comb = sorted([ptypes[a], ptypes[b]])
                if type_comb == sorted(self.types):
                    p_a = positions[a]
                    p_b = positions[b]
                    v_a = velocities[a]
                    v_b = velocities[b]

                    # world's most complicated way to get particle masses...
                    m_a = data.particles.particle_types.types[ptypes[a] - 1].mass
                    m_b = data.particles.particle_types.types[ptypes[b] - 1].mass

                    # calculate positional data
                    com, d_ab, d_ba, R_a, R_b = calc_positional(m_a, m_b,
                                                               p_a, p_b)

                    data.particles_['d'][a] = d_ab
                    data.particles_['d'][b] = d_ba
                    data.particles_['R'][a] = R_a
                    data.particles_['R'][b] = R_b

                    # calculate the velocities
                    v_com, v_vib_a, v_vib_b, v_rot_a, v_rot_b = calc_velocities(
                        m_a, m_b, d_ab, d_ba, v_a, v_b)

                    data.particles_['v_com'][a] = v_com
                    data.particles_['v_com'][b] = v_com

                    data.particles_['v_vib'][a] = v_vib_a
                    data.particles_['v_vib'][b] = v_vib_b
                    
                    data.particles_['v_rot'][a] = v_rot_a
                    data.particles_['v_rot'][b] = v_rot_b
        
        return 0


def calc_positional(m_a, m_b, p_a, p_b):
    """
    Compute positional data, i.e. center-of-mass, bond vector
    and bond radius vector.
    """
    # center-of-mass
    com = (m_a * p_a + m_b + p_b) / (m_a + m_b)
    # bond vector
    d_ab = p_b - p_a
    d_ba = p_a - p_b
    # bond radius vector
    R_a = p_a - com
    R_b = p_b - com

    return com, d_ab, d_ba, R_a, R_b 

def calc_velocities(m_a, m_b, d_ab, d_ba, v_a, v_b):
    """
    Compute the mode velocity vectors
    """
    ### center-of-mass
    v_com = (m_a * v_a + m_b * v_b) / (m_a + m_b)

    ### vibrational
    # bond unit vector
    dhat_ab = d_ab / norm(d_ab)
    dhat_ba = d_ba / norm(d_ba)
    # vibr vel
    v_vib_a = np.dot(v_a - v_com, dhat_ab) * dhat_ab
    v_vib_b = np.dot(v_b - v_com, dhat_ba) * dhat_ba

    ### rotational 
    v_rot_a = v_a - v_com - v_vib_a
    v_rot_b = v_b - v_com - v_vib_b

    return v_com, v_vib_a, v_vib_b, v_rot_a, v_rot_b







