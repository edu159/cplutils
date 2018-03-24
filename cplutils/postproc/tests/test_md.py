from cplutils.postproc.md import correlation
from cplutils.readers.lammps import chunk_reader
import pytest

class TestPostprocMd:
    # def test_correlation_atomistic(self, datadir):
    #     frames = chunk_reader(datadir[""], "molecule", order="frame")

    def test_correlation_mol(self, datadir):
        frames, time_idx = chunk_reader(str(datadir["molecule.corr"]), "molecule", order="frame")
        corr_func = correlation(frames)
        # 10 origins used
        assert len(corr_func) == 40

    def test_correlation_mol_origerror(self, datadir):
        frames, time_idx = chunk_reader(str(datadir["molecule.corr"]), "molecule", order="frame")
        with pytest.raises(Exception):
            corr_func = correlation(frames, orig_points=50)
        with pytest.raises(Exception):
            corr_func = correlation(frames, orig_points=100)
