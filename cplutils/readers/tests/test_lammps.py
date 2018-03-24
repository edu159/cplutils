from cplutils.readers.lammps import log_reader, chunk_reader
import pytest

# content of test_class.p
class TestLogReader:
    def test_1section(self, datadir):
        output = log_reader(str(datadir["log_1section.lammps"]),\
                            ["EQUILIBRATION"])
        section = output["EQUILIBRATION"]
        assert set(section["fields"].keys()) ==\
               set(["v_control_slab_density", "c_system_temp",\
                    "v_control_slab_press", "v_ctrl_slab_no_particles",\
                    "v_fluid_press", "v_fluid_width"]) 
        assert len(section["time"]) == 3

    def test_1section_nosections(self, datadir):
        output = log_reader(str(datadir["log_1section.lammps"]))
        assert set(output.keys()) == set(["run0", "run1"])
        assert len(output["run0"]["fields"]) == 5
        assert len(output["run0"]["fields"]["Temp"]) == 2 
        assert len(output["run1"]["fields"]) == 6
        assert len(output["run1"]["fields"]["v_ctrl_slab_no_particles"]) == 3

    def test_2sections(self, datadir):
        output = log_reader(str(datadir["log_2sections.lammps"]), ["EQUILIBRATION", "COUPLED"])
        assert set(output.keys()) == set(["EQUILIBRATION", "COUPLED"])
        assert len(output["EQUILIBRATION"]["fields"]) == 6
        assert len(output["EQUILIBRATION"]["fields"]["v_ctrl_slab_no_particles"]) == 2001
        assert len(output["COUPLED"]["fields"]) == 9
        assert len(output["COUPLED"]["fields"]["v_wall_velcom"]) == 10001


    def test_2sections_nosections(self, datadir):
        output = log_reader(str(datadir["log_2sections.lammps"]))
        assert set(output.keys()) == set(["run0", "run1"])
        assert len(output["run0"]["fields"]) == 6
        assert len(output["run0"]["fields"]["v_ctrl_slab_no_particles"]) == 2001
        assert len(output["run1"]["fields"]) == 9
        assert len(output["run1"]["fields"]["v_wall_velcom"]) == 10001

    def test_badrow(self, datadir):
        with pytest.raises(Exception):
            output = log_reader(str(datadir["log_badrow.lammps"]), ["COUPLED"])

    def test_incomplete(self, datadir):
        output = log_reader(str(datadir["log_incomplete.lammps"]), ["COUPLED"])
        assert len(output["COUPLED"]["fields"]) == 9
        assert len(output["COUPLED"]["fields"]["v_wall_velcom"]) == 21

    def test_1run(self, datadir):
        output = log_reader(str(datadir["log_1run.lammps"]), ["COUPLED"])


class TestChunkReader:
    def test_molecule_chunk_frame(self, datadir):
        output = chunk_reader(str(datadir["molecule.chunk"]), "molecule", order="frame")

    def test_molecule_chunk_row(self, datadir):
        output = chunk_reader(str(datadir["molecule.chunk"]), "molecule", order="row")

    def test_trajectory_frame(self, datadir):
        output = chunk_reader(str(datadir["trajectory-xyz.chunk"]), "trajectory-xyz", order="frame")

    def test_trajectory_row(self, datadir):
        output = chunk_reader(str(datadir["trajectory-xyz.chunk"]), "trajectory-xyz", order="row")

    def test_spatial1D(self, datadir):
        output = chunk_reader(str(datadir["spatial-1d.chunk"]), "spatial-1d", L=[10])

    def test_spatial3D(self, datadir):
        output = chunk_reader(str(datadir["spatial-3d.chunk"]), "spatial-3d", L=[10,10,10])

