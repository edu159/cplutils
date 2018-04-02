from cplutils.readers.lammps import log_reader, chunk_reader
import pytest
import numpy as np

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
        output = chunk_reader(str(datadir["spatial-1d.chunk"]), "spatial-1d")
        assert len(output["axis"]) == 1
        assert list(output["axis"][0]) == list(np.linspace(1.0, 59.0, 30))

    def test_spatial3D(self, datadir):
        output = chunk_reader(str(datadir["spatial-3d.chunk"]), "spatial-3d")
        assert len(output["axis"]) == 3
        assert np.isclose(output["axis"][0], np.linspace(3.61, 32.49, 5)).all()
        assert np.isclose(output["axis"][1], np.linspace(4.79272, 62.9056, 11)).all()
        assert np.isclose(output["axis"][2], np.linspace(3.61, 32.49, 5)).all()

    def test_ave_time_array_row(self, datadir):
        output = chunk_reader(str(datadir["ave-time-array.chunk"]), "ave/time-array", order="row")

    def test_ave_time_array_frame(self, datadir):
        output = chunk_reader(str(datadir["ave-time-array.chunk"]), "ave/time-array", order="frame")
