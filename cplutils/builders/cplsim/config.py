#from __future__ import with_statement
import yaml
import os
import re


class ConfigFile:
    def __init__(self):
        self.ALLOWED_SECTIONS = {"CFD-OPTIONS": OpenFOAMConfigSection,
                                 "COUPLING-OPTIONS": CPLConfigSection}
                               #  "MOLECULAR-OPTIONS": "", 
                               #  "COUPLING-OPTIONS": ""}
        self.config_fname = ""
        self.loaded = False
        self.config_data = {}
        self.sections ={}

    def load(self, fname="config.cfg"):
        try:
            with open(fname, 'r') as configfile:
                self.config_data = yaml.load(configfile)
        except Exception as error:
            print "Error in configuration file:"
            print "\t", error
        self._load_sections()
        self.config_fname = fname
        self.loaded = True

    def _load_sections(self):
        for section_name, section_opts in  self.config_data.items():
            try:
                section_class =  self.ALLOWED_SECTIONS[section_name]
                self.sections[section_name] = section_class(section_name, section_opts)
            except Exception as error:
                print "Error: section not found ", error 

            
    def add_section(self, section, opts):
        self.config_data[section] = opts

class ConfigSection:
    def __init__(self, section_name, section_opts, required_opts):
        self.section_name = section_name
        self.section_opts = section_opts
        self.required_opts = required_opts
        self._check_opts()

    def _check_opts(self):
        def _check_opts_recursive(opts_dict, opts_types_dict):
            for opt_name, opt_val in opts_dict.items():
                try:
                    types = opts_types_dict[opt_name]
                except Exception as error:
                    print "Option not found: ", error
                if type(types) is dict:
                    _check_opts_recursive(opts_dict[opt_name], opts_types_dict[opt_name])
                elif type(types) is tuple:
                    if type(opt_val) not in types:
                        print "Type not correct in option ", opt_name, " value: ", opt_val
                elif type(types) is str and opt_val not in types.split('|'):
                    print "Type incorrect in option ", opt_name, " value: ", opt_val
        _check_opts_recursive(self.section_opts, self.required_opts)

class OpenFOAMConfigSection(ConfigSection):
    def __init__(self, section_name, section_opts):
        required_opts = {"engine": (str,), 
         "version": (str,), 
         "solver": (str,), 
         "mesh": {"ncx": (int,), "ncy": (int,), "ncz": (int,)},
         "domain" : {"lx": (int, float), "ly": (int, float), "lz": (int, float)}, 
         "time" : {"start": (float, int), "end": (float, int), "delta": (float, int), "write-interval": (float, int)}, 
         "parallel": {"npx": (int), "npy": (int), "npz": (int)},
         "liquid-properties": {"density": (float, str), "viscosity": (float, str)},
         "boundary-field": {"velocity": {"mode": "fixed-value|zero-gradient|constant-gradient", 
             "x": (float, int, str),
             "y": (float, int, str), 
             "z": (float, int, str)},
             "pressure": {"mode": "fixed-value|zero-gradient|constant-gradient",
                 "x": (float, int, str),
                 "y": (float, int, str), 
                 "z": (float, int, str)}},
             "coupling-opts": {"stress-compute": {"mode": "surface|cell-centre"}}}
        ConfigSection.__init__(self, section_name, section_opts, required_opts)

class CPLConfigSection(ConfigSection):
    def __init__(self, section_name, section_opts):
        required_opts = {"timestep-ratio": (int), 
                "units": "lj|real", 
                "coupling-unit": {"x-unit": (float), "y-unit": (float), "z-unit": (float)}, 
                "const-region": {"lower-cell": (int), "upper-cell": (int)},
                "boundary-region": {"lower-cell": (int), "upper-cell": (int)}, 
                "coupling-region": {"cells": (int)}}
        ConfigSection.__init__(self, section_name, section_opts, required_opts)
 



class ConfigGenerator:
    def __init__(self, config_file, section_name):
        self.config_file = config_file
        self.template_base_dir = "templates"
        self.template_dir = ""
        self.files = []
        self.section = self.config_file.sections[section_name]
        self.computed_opts = {}
        self.computed_opts["computed"] = {}
        self._compute_options()

    def generate(self, output_dir):
        for f in self.files:
            try:
                with open(os.path.join(self.template_dir, f), 'r') as opts_file:
                    self._replace_opts(opts_file)
            except Exception as error:
                print error


    def _replace_opts(self, file):
        lines = file.readlines()                                                                                                                                                                                                 
        for line_idx in xrange(len(lines)):
            line_opts = re.findall(r'\$\[([a-zA-Z0-9\-\.]+?)\]', lines[line_idx])
            for opt in line_opts:
                opt_val = self._dots2dict(opt)
                lines[line_idx] = lines[line_idx].replace("$[" + opt + "]", str(opt_val))
        print "".join(lines)

    def _dots2dict(self, opt):
        keys = opt.split('.')
        opt_val = None
        opt_section_val = self.section.section_opts
        opt_computed_val = self.computed_opts 
        opt_type = self.section.required_opts
        try:
            for key in keys:
                opt_type = opt_type[key]
                opt_section_val = opt_section_val[key]
            opt_val = self._postproc_opt(opt_section_val, opt_type)
            #print opt_section_val,":", opt_type
        except Exception as error:
            try:
                for key in keys:
                    opt_computed_val = opt_computed_val[key]
                opt_val = opt_computed_val
            except Exception as error:
                print "Error: ", error
        #print opt, ":", opt_val
        return opt_val

    def _postproc_opt(self, opt_in, opt_type):
        return opt_in

    def _compute_options():
        pass



#TODO: Replace str for func type which parse the str

class OpenFOAMConfigGenerator(ConfigGenerator):
    def __init__(self, config_file):
        ConfigGenerator.__init__(self, config_file, "CFD-OPTIONS")
        self.foam_version = self.section.section_opts["version"]
        self.template_dir = os.path.join(self.template_base_dir, "openfoam/" + self.foam_version)
        self.files = ["0/p", "0/U", "constant/transportProperties", "constant/polyMesh/blockMeshDict",
                "system/controlDict", "system/decomposeParDict"]

    def _compute_options(self):
        # Compute the subdomains
        npx = self.section.section_opts["parallel"]["npx"] 
        npy = self.section.section_opts["parallel"]["npy"] 
        npz = self.section.section_opts["parallel"]["npz"] 
        self.computed_opts["computed"]["num-subdomains"] = npx * npy * npz
        # Compute kinematic viscosity
        density = self.section.section_opts["liquid-properties"]["density"] 
        viscosity = self.section.section_opts["liquid-properties"]["viscosity"] 
        self.computed_opts["computed"]["nu"] = viscosity/density

    def _postproc_opt(self, opt_in, opt_type):
        opt_out = opt_in
        # If the option is a choice one then use camel case as OpenFOAM does
        if type(opt_type) is str:
            opt_out = [word.title() for word in opt_in.split('-')]
            opt_out[0] = opt_out[0].lower()
            opt_out = "".join(opt_out)
        return opt_out


class CPLConfigGenerator(ConfigGenerator):
    def __init__(self, config_file):
        ConfigGenerator.__init__(self, config_file, "COUPLING-OPTIONS")
        self.template_dir = os.path.join(self.template_base_dir, "cpl")
        self.files = ["COUPLER.in"]

    def _compute_options(self):
        cfd_section = self.config_file.sections["CFD-OPTIONS"]
        self.computed_opts["computed"]["density-cfd"] = cfd_section.section_opts["liquid-properties"]["density"]
        self.computed_opts["computed"]["icmin-olap"] = 1
        self.computed_opts["computed"]["icmax-olap"] = cfd_section.section_opts["mesh"]["ncx"]
        self.computed_opts["computed"]["jcmin-olap"] = 1
        self.computed_opts["computed"]["jcmax-olap"] = self.section.section_opts["coupling-region"]["cells"]
        self.computed_opts["computed"]["kcmin-olap"] = 1
        self.computed_opts["computed"]["kcmax-olap"] = cfd_section.section_opts["mesh"]["ncz"]
        self.computed_opts["computed"]["icmin-cnst"] = 1
        self.computed_opts["computed"]["icmax-cnst"] = cfd_section.section_opts["mesh"]["ncx"]
        self.computed_opts["computed"]["jcmin-cnst"] = self.section.section_opts["const-region"]["lower-cell"]
        self.computed_opts["computed"]["jcmax-cnst"] = self.section.section_opts["const-region"]["upper-cell"]
        self.computed_opts["computed"]["kcmin-cnst"] = 1
        self.computed_opts["computed"]["kcmax-cnst"] = cfd_section.section_opts["mesh"]["ncz"]
        


#TODO: Add temperature and options based on the solver, eg. if density or temperature options has to be
       # taken into account

class CaseBuilder:
    def __init__(self, path, case_name):
        self.case_name = case_name

    def _build(self):
        pass


if __name__ == "__main__":
    case_folder = "/home/edu/Desktop/repositories/cplsim/case"
    config_file = ConfigFile()
    config_file.load()
    c = OpenFOAMConfigGenerator(config_file)
    c.generate(case_folder)
    d = CPLConfigGenerator(config_file)
    d.generate(case_folder)

