from cplutils.readers.openfoam import field_reader
import pytest

class TestFoamReader:
    def test_field_reader_boundary(self, datadir):
        boundary_patches = boundary_patches = [("CPLReceiveMD",(0, -1, 0)),\
                                               ("movingWall",(0, 1, 0))]
        output = field_reader(str(datadir["openfoam_case"]), 
                              "U", boundary=True,
                              boundary_patches=boundary_patches)

    def test_field_reader_noboundary(self, datadir):
        output = field_reader(str(datadir["openfoam_case"]), "U")
