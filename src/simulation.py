import os

def find(path, filename):
    """
    Find first matching file in a directory according to filename.
    
    Returns:
    relative path to the first match of the file searching.
    """
    for root, dirs, files in os.walk(path):
        if filename in files:
            return os.path.join(root, filename)


class Simulation:
    """
    This defines a simulation class, which holds all relevant 
    files of a LAMMPS simulation. Child classes for all the analysis 
    methods will be created, which then inherent the properties of 
    this class.
    All of them are optional to state, but an error will be raised 
    in the respective method if a property is missing that is needed to 
    perform the desired analysis.
    """

    def __init__(self, 
                 log_file = '', 
                 traj_file = '',
                 bond_file = '',
                 vacf_file = '', 
                 dipole_file = ''):
        
        """
        Constructor for class Simulation.
        Uses the find method to find the path to all the properties
        defined for this class.
        
        Parameters:
        log_file (str):     log file
        traj_file (str):    trajectory / dump file
        bond_file (str):    reaxff bond file
        vacf_file (str):    file containing the velocity autocorrelation 
                            function as a function of time
        dipole_file (str):  file containing the dipole moment 
                            function as a function of time
        plot (bool):        determines if a plot from the data is created 
                            or not (Default: False)
        """
        self.log_file = find('.', log_file)
        self.traj_file = find('.', traj_file)
        self.bond_file = find('.', bond_file)
        self.vacf_file = find('.', vacf_file)
        self.dipole_file = find('.', dipole_file)


