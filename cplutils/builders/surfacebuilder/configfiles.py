from ase.io import read
import yaml
import shutil
import os

class SurfaceSection:
    pass

class PotentialSection:
    pass

class SurfacePotentialFile:
    def __init__(self, file_path):
        self.ALLOWED_SECTIONS = {"SURFACE": SurfaceSection, 
                                 "POTENTIAL": PotentialSection}
        self.path = file_path
        self.loaded = False
        self.data = {}
        self.sections ={}

    def load(self):
        try:
            with open(self.path, 'r') as surfacefile:
                self.data = yaml.safe_load(surfacefile)
        except IOError as e:
            raise Exception("Problem opening '%s' potential file - %s." % (self.path, e.strerror))
        except Exception as error:
            raise Exception("Parsing error in potential file '%s': \n" % self.path + str(error))
        self._load_sections()
        self.loaded = True


    def _load_sections(self):
        pass
            
    def __getitem__(self, key):
        if self.loaded:
            return self.data[key]
        else:
            raise Exception("Potential file '%s' not loaded." % self.path)


class SurfaceDB:
    def __init__(self, path=None):
        self.infofname = "surfaces.yaml"
        self.dirname = "surfaces-db"
        if path is None:
            self.path = os.getenv("HOME")
        else:
            self.path = path
        self.path = os.path.join(os.path.abspath(self.path), self.dirname)
        self.path_data = os.path.join(self.path, "data")
        self.loaded = False
        self.surfaces = None

    @staticmethod
    def create_db(path):
        db_path = os.path.join(path, "surfaces-db")
        try:
            os.mkdir(db_path)
            os.mkdir(os.path.join(db_path, "data"))
            open(os.path.join(db_path, "surfaces.yaml"), 'w').close()
        except Exception as err:
            raise Exception("Error creating database - {}.".format(str(err)))


    def load(self):
        try:
            with open(os.path.join(self.path, self.infofname), 'r') as infofile:
                self.surfaces = yaml.safe_load(infofile)
        except IOError as e:
            raise Exception("Problem opening '%s' file - %s." % (self.infofname, e.strerror))
        except Exception as error:
            raise Exception("Parsing error in '%s': \n" % self.infofname + str(error))
        self.loaded = True

    def _copy_all_files(self, path, dest):
        for file_name in os.listdir(path):
            fpath = os.path.join(path, file_name)
            shutil.copy(fpath, dest)
 

    def copy_surfacedata(self, name, topology, potential, dest): 
        data_dir = os.path.join(self.path, "data")
        surface_dir = os.path.join(data_dir, name)
        topology_dir = os.path.join(surface_dir, "topologies")
        potential_dir = os.path.join(surface_dir, "potentials")
        # The two files to copy
        topology_fname = topology + ".xyz"
        potential_fname = "{}.yaml".format(potential)
        shutil.copy(os.path.join(topology_dir, topology_fname), dest)
        shutil.copy(os.path.join(potential_dir, potential_fname), dest)
        
    def _add_to_existent(self, surface_name, topology_name, atom_types, dimensions, nof_atoms, potential_fname):
        surface = self.surfaces[surface_name]
        assert atom_types == surface["atom-types"]
        assert topology_name not in surface["topologies"].keys()
        surface["topologies"][topology_name] = {}
        surface["topologies"][topology_name]["dimensions"] = dimensions
        surface["topologies"][topology_name]["nof-atoms"] = nof_atoms
        if potential_fname is not None:
            if potential_fname in surface["potentials"]:
                raise Exception("Error: Potential '{}' already exists for surface '{}'.".format(potential_fname, surface_name))
            surface["potentials"].append(potential_fname)

    def _add_new(self, surface_name, topology_name, atom_types, dimensions, nof_atoms, potential_fname):
        if not self.surfaces:
            self.surfaces = {}
        if potential_fname is None:
            raise Exception("Error: Potential file not supplied for a new surface. Plese specify one.")
        self.surfaces[surface_name] = {}
        surface = self.surfaces[surface_name]
        surface["atom-types"] = atom_types
        surface["topologies"] = {}
        surface["topologies"][topology_name] = {}
        surface["topologies"][topology_name]["dimensions"] = dimensions
        surface["topologies"][topology_name]["nof-atoms"] = nof_atoms
        surface["potentials"] = [potential_fname]



    def add(self, surface_name, topology_name, topology_file, potential_file=None):
        assert self.loaded
        topology = read(topology_file)
        dimensions = [float(topology.cell[0,0]), float(topology.cell[1,1]), float(topology.cell[2,2])]
        nof_atoms = topology.get_number_of_atoms() 
        atom_types = list(set(topology.get_chemical_symbols()))
        surface_dir = os.path.join(self.path_data, surface_name)
        topology_dir = os.path.join(surface_dir, "topologies")
        potential_dir = os.path.join(surface_dir, "potentials")
        if potential_file is not None:
            pf = SurfacePotentialFile(potential_file)
            pf.load()
            potential_fname = pf["POTENTIAL"]["name"]
        else:
            potential_fname = None
        if surface_name in self.surfaces.keys():
            self._add_to_existent(surface_name, topology_name, atom_types, dimensions, nof_atoms, potential_fname)
        else:
            self._add_new(surface_name, topology_name, atom_types, dimensions, nof_atoms, potential_fname)
            os.mkdir(surface_dir)
            os.mkdir(topology_dir)
            os.mkdir(potential_dir)
        # Copy topology and potential
        shutil.copy(topology_file, os.path.join(topology_dir, topology_name + ".xyz"))
        if potential_file is not None:
            shutil.copy(potential_file, os.path.join(potential_dir, potential_fname + ".yaml"))

    def save(self):
        with open(os.path.join(self.path, self.infofname), 'w') as infofile:
            yaml.dump(self.surfaces, infofile, default_flow_style=False)
    
    def delete_surfaces(self, mol_names):
        pass

    def imprt(self):
        pass

    def export(self):
        pass


