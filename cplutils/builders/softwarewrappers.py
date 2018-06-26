from subprocess import Popen, PIPE
import subprocess
import select
import time
import re
import os
import sys

class MTSystemWriter:
    def __init__(self, box_dims, molecules, system_name='system.lt'):
        self.box_xlo = box_dims[0]
        self.box_xhi = box_dims[3]
        self.box_ylo = box_dims[1]
        self.box_yhi = box_dims[4]
        self.box_zlo = box_dims[2]
        self.box_zhi = box_dims[5]
        self.molecules = molecules
        self.fname = system_name
        self.lines = []



    def _write_header(self):
        for mol in self.molecules:
            self.lines.append('import ' + mol[0] + '.lt' + '\n')

    def _write_molecules(self):
        m = 0
        for mol in self.molecules:
            self.lines.append('mol' + str(m) + ' = new ' + mol[0] + '[%d]' % mol[1] + '\n')
            m += 1

    def _write_box(self):
        self.lines.append('write_once("Data Boundary") {' + '\n')
        self.lines.append('    ' + str(self.box_xlo) + ' ' + str(self.box_xhi) + ' xlo xhi' + '\n')
        self.lines.append('    ' + str(self.box_ylo) + ' ' + str(self.box_yhi) + ' ylo yhi' + '\n')
        self.lines.append('    ' + str(self.box_zlo) + ' ' + str(self.box_zhi) + ' zlo zhi' + '\n')
        self.lines.append('}' + '\n')


    def build_system_file(self):
        self._write_header()
        self.lines.append('\n')
        self._write_molecules()
        self.lines.append('\n')
        self._write_box()

    def save(self, path='.'):
        with open(os.path.join(path, self.fname), 'w') as sysmolfile:
            sysmolfile.writelines(self.lines)



class Moltemplate():
    def __init__(self, ftype="xyz", atomstyle="full", topology_fname=None, 
                 template_fname="system.lt", exec_path="moltemplate.sh"):
        if topology_fname is None:
            self.topology_fname = "system.%s" % ftype
        else:
            self.topology_fname = topology_fname
        self.ftype = ftype
        self.atomstyle = atomstyle
        self.template_fname = template_fname
        self.exec_path = exec_path

    def run(self):
        output = ""
        success = False
        try:
            output = subprocess.check_output([self.exec_path, "-%s" % self.ftype, 
                                    self.topology_fname, "-atomstyle", "full",
                                    self.template_fname], stderr=subprocess.STDOUT)
            success = True
        except subprocess.CalledProcessError as error:
            #TODO: Delete junk intermediate files from moltemplate
            output = error.output
            raise Exception("Error:\n" + str(error.output))
        finally:
            with open("moltemplate.log", "w") as log:
                log.write("\n ---- moltemplate.sh log ---\n\n")
                log.writelines(output)

        # Clean unused atomtypes
        # TODO: Refactor into a function/class to run and log outputs (packmol, moltemplate)
        if success:
            try:
                output = subprocess.check_output(["cleanup_moltemplate.sh"], stderr=subprocess.STDOUT)
            except subprocess.CalledProcessError as error:
                output = error.output
                raise Exception("Error:\n" + str(error.output))
            finally:
                with open("moltemplate.log", "a") as log:
                    log.write("\n ---- cleanup_moltemplate.sh log ---\n\n")
                    log.write(output)


class Packmol():
    def __init__(self, packmol_fname="system.packmol", exec_path="packmol"):
        self.exec_path = exec_path
        self.packmol_fname= packmol_fname
        self.max_time = 15 # In minutes
        self.current_const_violation = None
        self.current_dist_violation = None
        self.box_tol = 2.0
        self.output_file = self._get_output_fname()

    def _constrain_violation(self, line):
        if "Maximum violation of target distance:" in line:
            print line
            print line.split(" ")[-1]
            self.current_dist_violation = float(line.split(" ")[-1])
            print "Distance violation: ", self.current_dist_violation
        elif "Maximum violation of the constraints:" in line:
            self.current_const_violation = float(line.split(" ")[-1])
            print "Constrain violation: ", self.current_const_violation
        elif "ERROR:" in line:
            raise Exception("Error in Packmol. Check packmol.out for details.")

        if self.current_const_violation is not None:
            if self.current_const_violation < self.box_tol:
                return False
        return True

    def _get_output_fname(self):
        with open(self.packmol_fname, 'r') as packmol_fname:
            lines = packmol_fname.readlines()
            for l in lines:
                if re.search("^[ ]*outpu",l) is not None:
                    return [c for c in l.split(" ") if c][-1][:-1]


    def _system_outputfile_exists(self):
        return os.path.exists(self.output_file)

    def _wait_for_output(self):
        fsize_ini = os.path.getsize(self.output_file)
        while True:
            time.sleep(1)
            fsize_curr = os.path.getsize(self.output_file)
            if fsize_ini == fsize_curr:
                break
            else:
                fsize_ini = fsize_curr

    def run(self):
        with open(self.packmol_fname, 'r') as packmol_fname:
            p = Popen([self.exec_path], stdin=packmol_fname, 
                      stdout=PIPE, stderr=PIPE)
            s=select.poll()
            s.register(p.stdout,select.POLLIN)
            time_ini = time.time()
            running = True
            while running:
                if (abs(time_ini - time.time()) / 60.0) < self.max_time: 
                    if s.poll(1):
                        line = p.stdout.readline()
                        if not line:
                            p.poll()
                            if p.returncode is not None:
                                running = False
                        else:
                            try:
                                if not self._constrain_violation(line) and\
                                       self._system_outputfile_exists():
                                    self._wait_for_output()
                                    running = False
                                    p.kill()
                            except Exception as error:
                                p.kill()
                                running = False
                                raise 

class PackmolSystemWriter:
    def __init__(self, box_dims, molecules, packmol_tol=1.5, seed=-1, ftype="xyz", system_name='system'):
        self.molecules = molecules
        self.box_dims = box_dims
        self.system_name = "%s.%s" % (system_name, ftype)
        self.packmol_fname ="%s.%s" % (system_name, "packmol") 
        self.lines = []
        self.packmol_tol = packmol_tol
        self.box_tol = 2.0
        self.ftype = ftype
        self.seed = seed

    def _apply_box_tol(self, box_dims):
        box_dims_out = list(box_dims)
        #TODO: Check box_dims are inside self.box_dims
        for i in xrange(3):
            delta = abs(box_dims[i] - self.box_dims[i])
            if delta < self.box_tol:
                print delta, self.box_tol
                box_dims_out[i] += self.box_tol-delta

        for i in xrange(3,6):
            delta = abs(box_dims[i] - self.box_dims[i])

            if delta < self.box_tol:
                print delta, self.box_tol
                box_dims_out[i] -= self.box_tol-delta

        return box_dims_out

    def _write_header(self):
        #TODO: Write comments here
        self.lines.append('tolerance %s\n' % self.packmol_tol)
        self.lines.append('writeout %d\n' % 1)
        # self.lines.append('precision %s\n' % self.box_tol)
        # Necessary as packmol has a default limit on side lenghts
        if max(self.box_dims) > 1000:
            self.lines.append('sidemax %s\n' % max(self.box_dims) + self.box_tol)
        self.lines.append('filetype %s\n' % self.ftype)
        self.lines.append('output %s\n' % self.system_name)
        self.lines.append('seed %s\n' % self.seed)

    def _write_molecule(self, molecule):
        self.lines.append('structure %s.%s\n' % (molecule["name"], self.ftype))
        self.lines.append('    number %s\n' % molecule["number"])

        if molecule["const_type"] == "box":
            if molecule["const"] is None:
                box_dims = self._apply_box_tol(self.box_dims)
            else:
                box_dims = tuple(self._apply_box_tol(molecule["const"]))
            self.lines.append('    inside box %s %s %s %s %s %s\n' % box_dims)
        elif molecule["const_type"] == "fixed":
                self.lines.append('    fixed %s %s %s %s %s %s\n' %\
                                  tuple(molecule["const"]))
        else:
            raise Exception("Constrain type in % molecule not recognised." % molecule["name"])
        self.lines.append('end structure')

    def build_system_file(self):
        self._write_header()
        for molecule in self.molecules:
            self.lines.append('\n')
            self._write_molecule(molecule)

    def save(self, path='.'):
        with open(os.path.join(path, self.packmol_fname), 'w') as sysfile:
            sysfile.writelines(self.lines)




