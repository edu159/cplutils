from ase.io import read
import os
import yaml
import numpy as np
from ase.neighborlist import neighbor_list
from configfiles import SurfacePotentialFile

# print a.get_atomic_numbers()
# print a.get_chemical_symbols()
# print a.get_positions()
#
class MTLennardJones:
    def __init__(self, pot_params):
        self.atom_list = pot_params["atoms"]
        self.type = pot_params["type"]

    def get_atom_labels(self):
        return self.atom_list.keys()

    def __getitem__(self, key):
        return self.atom_list[key]

    def get_coeff_strings(self):
        strings = []
        s = ""
        for atom, params in self.atom_list.items():
            if self.type == "lj/cut":
                s = "pair_coeff @atom:%s @atom:%s lj/cut %f %f %f" %\
                    (atom, atom, params["epsilon"], params["sigma"], params["cutoff"])
                strings.append(s)
            else:
                raise Exception("Lennard Jones pair type no recongnised.")
        return strings

class MTHarmonic:
    def __init__(self, pot_params):
        self.bonds = pot_params

    def __getitem__(self, key):
        return self.atom_list[key]

    def get_coeff_strings(self):
        strings = []
        s = ""
        for bond, params in self.bonds.items():
            s = "bond_coeff @bond:%s harmonic %f %f" %\
                    (bond, params["kbond"], params["rbond"])
            strings.append(s)
        return strings



class MTSurfaceWriter:
    def __init__(self, surface_name, surface_path, potential_path, fname="surface.lt"):
        # Allowed potentials and mapper
        spf = SurfacePotentialFile(potential_path)
        spf.load()
        self.potentials = spf["POTENTIAL"]
        self.atoms = read(surface_path)
        self.ALLOWED_PAIR = ["LennardJones"]
        self.ALLOWED_BOND = ["Harmonic"]
        self.POTENTIALS = {"LennardJones": MTLennardJones,
                           "Harmonic": MTHarmonic,
                          }
        # All the elements present in the surface
        self.atom_symbols = list(set(self.atoms.get_chemical_symbols()))
        self.atom_masses = list(set(self.atoms.get_masses()))
        self.atom_charges = [self.potentials["Coulombic"][a] for a in self.atom_symbols]
        print self.atom_symbols, self.atom_masses, self.atom_charges
        # self.atom_charges = list(set(self.atoms.get_chemical_symbols()))
        # TODO:Maybe hybrid types can be implemented at some point
        self.pair_potentials = self._init_potentials(self.ALLOWED_PAIR)
        self.bond_potentials = self._init_potentials(self.ALLOWED_BOND)
        self.surface_name = surface_name
        self.lines = []
        self.fname = fname
        self.L = self.atoms.get_cell_lengths_and_angles()[0:3]
        self.atom_positionsall = self.atoms.get_positions()


    def _init_potentials(self, pot_allowed):
        result = []
        # Initialize and return a list of present potentials
        for potential, params in self.potentials.items():
            if potential in pot_allowed:
                self.POTENTIALS[potential] = self.POTENTIALS[potential](params)
                result.append(potential)
        return result


    def _write_header(self):
        self.lines.append("%s {" % self.surface_name)
        self.lines.append('\n')

    def _write_footer(self):
        self.lines.append("} #%s " % self.surface_name)

    def _write_charges(self):
        self.lines.append('  write("In Charges") {\n')
        for i, c in enumerate(self.atom_charges):
            s = "    set type @atom:%s charge %f\n" % (self.atom_symbols[i], c)
            self.lines.append(s)
        self.lines.append('  }\n')
        self.lines.append('\n')

    def _write_masses(self):
        self.lines.append('  write_once("Data Masses") {\n')
        for i, m in enumerate(self.atom_masses):
            s = "    @atom:%s %f\n" % (self.atom_symbols[i], m)
            self.lines.append(s)
        self.lines.append('  }\n')
        self.lines.append('\n')

    def _write_coeffs(self):
        self.lines.append('  write("In Settings") {\n')
        self.lines.append('    # Surface potential coefficients\n')
        lines = []
        present_potentials = self.pair_potentials + self.bond_potentials
        for potential in present_potentials:
            pot_params = self.potentials[potential]
            pot = self.POTENTIALS[potential]
            lines.extend(pot.get_coeff_strings())
        lines = ['    ' + l + '\n' for l in lines]
        self.lines.extend(lines)
        self.lines.append('  }\n')
        self.lines.append('\n')

    def _write_data_atoms(self):
        self.lines.append('\n')
        self.lines.append('  # atomID   molID  atomTyle  charge  X   Y   Z\n')
        self.lines.append('  write("Data Atoms") {\n')
        for atom in self.atoms:
            atom_col = '    $atom:%d' % atom.index
            mol_col = '$mol:.'
            atomid_col = '@atom:%s' % atom.symbol
            charge_col = '0.0'
            atom_pos = atom.position
            x_col = str(atom_pos[0])
            y_col = str(atom_pos[1])
            z_col = str(atom_pos[2])
            row = '{0:<8} {1:<8} {2:<8} {3:<4} {4:<16} {5:<16} {6:<16} '\
                  .format(atom_col, mol_col, atomid_col,
                          charge_col, x_col, y_col, z_col)
            self.lines.append(row + '\n')
        self.lines.append('  }\n')
        self.lines.append('\n')

    def _distance(self, atom1, atom2, L, pbc=[False, False, False]):
        atom_pos = self.atom_positionsall
        delta = atom_pos[atom1] - atom_pos[atom2]
        pos_corrected = np.abs(delta)
        Li = 0.0
        for i in xrange(3):
            if pbc[i] and (abs(delta[i]) > L[i]/2.0):
                if delta[i] < 0.0:
                    Li = L[i]
                else:
                    Li = -L[i]
                pos_corrected[i] = abs(delta[i] +  Li)
        return pos_corrected[0]**2+pos_corrected[1]**2+pos_corrected[2]**2
        # return pos_corrected[0]+pos_corrected[1]+pos_corrected[2]

    def _generate_bonds(self, group1, group2, neighbors, dmin, dmax):
        ai = neighbors[0]
        aj = neighbors[1]
        distanceij = neighbors[2]
        bonds = set()
        atom_nof_neighs = np.bincount(ai)
        neigh_idx = np.cumsum(atom_nof_neighs)
        neigh_idx[1:] = neigh_idx[:-1]
        neigh_idx[0] = 0
        # print neigh_idx
        print "starts..."
        import time
        init = time.time()
        # If the groups are the same then things speed up
        groups_equal = (group1 == group2)
        for atom1 in group1:
            atom1_nof_neighs = atom_nof_neighs[atom1]
            if atom1_nof_neighs > 0:
                idx = neigh_idx[atom1]
                atom1_neighs = aj[idx:idx+atom1_nof_neighs]
                for idx2, atom2  in enumerate(atom1_neighs):
                    atom_in_group2 = True
                    if not groups_equal and group1atom2 not in group2:
                        atom_in_group2 = False
                    if atom_in_group2:
                        d = distanceij[idx+idx2]
                        if d >= dmin and d <= dmax: 
                            bond_key = tuple(sorted([atom1, atom2]))
                            bonds.add(bond_key)
        print time.time() - init
        return sorted(list(bonds))


    # Deprecated, very slow
    def _generate_bonds2(self, group1, group2, neighbors, distances):
        bonds = set()
        d = 0.0
        bond_key = ()
        for atom1 in group1:
            for atom2 in group2:
                if atom1 != atom2:
                    bond_key = tuple(sorted([atom1, atom2]))
                    if bond_key not in bonds:
                        L = self.L
                        pbc = self.atoms.pbc
                        atom_pos = self.atom_positionsall
                        delta = atom_pos[atom1] - atom_pos[atom2]
                        pos_corrected = np.abs(delta)
                        Li = 0.0
                        for i in xrange(3):
                            if pbc[i] and (abs(delta[i]) > L[i]/2.0):
                                if delta[i] < 0.0:
                                    Li = L[i]
                                else:
                                    Li = -L[i]
                                pos_corrected[i] = abs(delta[i] +  Li)

                        d = pos_corrected[0]**2+pos_corrected[1]**2+pos_corrected[2]**2
                        if d >= distances[0]**2 and d <= distances[1]**2:
                            bond_key = tuple(sorted([atom1, atom2]))
                            bonds.add(bond_key)
        return sorted(list(bonds))


    def _write_data_bonds(self, bonds):
        self.lines.append("  write('Data Bonds') {\n")
        for bond_name, bond_values in bonds.items():
            for bond_idx in bond_values:
                bond_col = '    $bond:%d-%d' % (bond_idx[0], bond_idx[1])
                bondname_col = '@bond:%s' % bond_name
                at1_col = '$atom:%d' % bond_idx[0]
                at2_col = '$atom:%d' % bond_idx[1]
                row = '{0:<12} {1:<10} {2:<12} {3:<12}'\
                      .format(bond_col, bondname_col, at1_col, at2_col)
                self.lines.append(row + '\n')
        self.lines.append('  }\n')
        self.lines.append('\n')

    def _get_all_bonds(self):
        bonds_all = {}
        group1 = []
        group2 = []
        for potential in self.bond_potentials:
            for bond, params in self.potentials[potential].items():
                if bond in bonds_all:
                    raise Exception("Error: Bond names have to be unique!")
                group1 = [a.index for a in self.atoms if a.symbol in params["creation"]["group1"]]
                group2 = [a.index for a in self.atoms if a.symbol in params["creation"]["group2"]]
                neighbors = neighbor_list("ijd", self.atoms, params["creation"]["distance"][1])
                distances = params["creation"]["distance"]
                bonds = self._generate_bonds(group1, group2, neighbors, distances[0], distances[1])
                bonds_all[bond] = bonds
        print bonds_all.keys(), len(bonds_all.values()[0]), len(bonds_all.values()[1])
        return bonds_all
        # print set(bonds_all["bond1"]).intersection(set(bonds_all["bond2"]))

    def build_surface_file(self):
        self._write_header()
        self._write_data_atoms()
        self._write_charges()
        self._write_masses()
        self._write_coeffs()
        bonds_all = self._get_all_bonds()
        self._write_data_bonds(bonds_all)
        self._write_footer()

    def save(self, path='.'):
        with open(os.path.join(path, self.fname), 'w') as surfacemolfile:
            surfacemolfile.writelines(self.lines)


