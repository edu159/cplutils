import os
import yaml
import shutil
class SystemSection:
    pass
class ParamsSection:
    pass

class MDSystemFile:
    def __init__(self, path='.', fname='mdsystem.yaml'):
        self.ALLOWED_SECTIONS = {"SYSTEM": SystemSection, 
                                 "PARAMETERS": ParamsSection}
        self.system_path = os.path.abspath(path)
        self.fname = fname
        self.path = os.path.join(self.system_path, fname)
        self.loaded = False
        self.params_data = {}
        self.sections ={}

    def load(self):
        try:
            with open(self.path, 'r') as mdsystemfile:
                self.params_data = yaml.load(mdsystemfile)
        except IOError as e:
            raise Exception("Problem opening '%s' file - %s." % (self.fname, e.strerror))
        except Exception as error:
            raise Exception("Parsing error in '%s': \n" % self.fname + str(error))
        self._load_sections()
        self.loaded = True


    def _load_sections(self):
        pass
            
    def __getitem__(self, key):
        if self.loaded:
            return self.params_data[key]
        else:
            raise Exception("File '%s' not loaded." % self.fname)


class MoleculeDB:
    def __init__(self, path=None):
        self.infofname = "molecules.yaml"
        self.dirname = "molecules-db"
        if path is None:
            self.path = os.getenv("HOME")
        else:
            self.path = path
        self.path = os.path.join(os.path.abspath(self.path), self.dirname)
        self.loaded = False

    def load(self):
        try:
            with open(os.path.join(self.path, self.infofname), 'r') as infofile:
                self.molecules = yaml.load(infofile)
        except IOError as e:
            raise Exception("Problem opening '%s' file - %s." % (self.infofname, e.strerror))
        except Exception as error:
            raise Exception("Parsing error in '%s': \n" % self.infofname + str(error))
        self.loaded = True

    def _copy_all_files(self, path, dest):
        for file_name in os.listdir(path):
            fpath = os.path.join(path, file_name)
            shutil.copy(fpath, dest)
 

    def copy_moldata(self, molname, potential, dest): 
        data_dir = os.path.join(self.path, "data")
        mol_dir = os.path.join(data_dir, "molecules")
        potential_dir = os.path.join(data_dir, "potentials")
        mol_dir = os.path.join(mol_dir, molname)
        moltemplate_dir = os.path.join(mol_dir, "moltemplate")
        # The three directories to copy
        topology_dir = os.path.join(mol_dir, "topology")
        potential_dir = os.path.join(potential_dir, potential)
        moltemplate_dir = os.path.join(moltemplate_dir, potential)
        self._copy_all_files(potential_dir, dest)
        self._copy_all_files(moltemplate_dir, dest)
        self._copy_all_files(topology_dir, dest)
        
    def save(self):
        pass
    
    def delete_molecules(self, mol_names):
        pass

    def impor(self):
        pass

    def export(self):
        pass


class SurfaceDB:
    def __init__(self, path=None):
        self.infofname = "surfaces.yaml"
        self.dirname = "surfaces-db"
        if path is None:
            self.path = os.getenv("HOME")
        else:
            self.path = path
        self.path = os.path.join(os.path.abspath(self.path), self.dirname)
        self.loaded = False

    def load(self):
        try:
            with open(os.path.join(self.path, self.infofname), 'r') as infofile:
                self.surfaces = yaml.load(infofile)
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
        moltemplate_dir = os.path.join(surface_dir, "moltemplate")
        # The two files to copy
        topology_fname = topology + ".xyz"
        moltemplate_fname = "%s_%s.lt" % (topology, potential)
        shutil.copy(os.path.join(topology_dir, topology_fname), dest)
        shutil.copy(os.path.join(moltemplate_dir, moltemplate_fname), dest)
        
    def save(self):
        pass
    
    def delete_surfaces(self, mol_names):
        pass

    def impor(self):
        pass

    def export(self):
        pass


