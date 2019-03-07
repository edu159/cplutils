import sys
import numpy as np
import glob
import os
try:
    from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile 
    from PyFoam.RunDictionary.ParsedBlockMeshDict import ParsedBlockMeshDict
except:
    print "Error: PyFoam package is required"
    sys.exit(1)

def timestep_foldername(step):
    if step >= 1e6:
        a = "{0:.16e}".format(step)
	print a.split('e')[0].rstrip('0').rstrip('.')+'e'+a.split('e')[1]
	return a.split('e')[0].rstrip('0').rstrip('.')+'e'+a.split('e')[1]
    else:
        return str(step)

def field_reader(casedir, field_name, time_steps=[], boundary=False,
                 boundary_patches=None):
    if boundary:
        if boundary_patches is None:
            raise Exception("Boundary patches have to be provided.")
    # Boundary fiel
    time_steps.sort()
    time_steps_all = [int(float(os.path.basename(d))) for d in\
                      glob.glob(os.path.join(casedir,"[1-9]*"))]
    time_steps_all.sort()
    if not time_steps:
        time_steps = time_steps_all
    else:
        for s in time_steps:
            if s not in time_steps_all:
                print "Timestep %s not present in openfoam case folder." % str(s)
                raise Exception()
    try:
        cell_cx = ParsedParameterFile(os.path.join(casedir, "0/ccx"))["internalField"]
        cell_cy = ParsedParameterFile(os.path.join(casedir, "0/ccy"))["internalField"]
        cell_cz = ParsedParameterFile(os.path.join(casedir, "0/ccz"))["internalField"]
    except Exception as e:
        pass
    try:
        blockMesh = ParsedBlockMeshDict(os.path.join(casedir,
                                        "constant/polyMesh/blockMeshDict"))
    except Exception as e:
        print e
        pass

    # Compute domain and mesh values
    ncx, ncy, ncz = tuple(blockMesh["blocks"][2])
    bounds = np.array(blockMesh.getBounds())
    Lx, Ly, Lz = bounds[1] - bounds[0]
    dx = Lx/ncx
    dy = Ly/ncy
    dz = Lz/ncz

    # Get field dimension
    time_dir = os.path.join(casedir, "%s/%s" % (str(0), field_name))
    field = ParsedParameterFile(time_dir)["internalField"]
    if type(field[0]) in [int, float]:
        field_size = 1
    else:
        field_size = len(field[0])
    #NOTE: Assumed in X and Z the boundaries are periodic
    #NOTE: Top boundary = movingWall, Low boundary = CPLReceiveMD
    #NOTE: Create an array assuming boundary=True, then slice at the end if False
    ncy_boundary = ncy + 2
    array = np.zeros((len(time_steps), ncx, ncy_boundary, ncz, field_size))
    
    def get_cell_internal(cell_coord, dx,dy,dz):
        return [int(cell_coord[0]/dx), int(cell_coord[1]/dy), int(cell_coord[2]/dz)]
    def get_cell_boundary():
        pass

    def patch2slice(patch):
        if patch == (None, None, None):
            return (slice(None), slice(None), slice(None))
        else:
            return [int(-(c+1) / 2.0) if c != 0 else slice(None) for c in patch]


    def field2array(field, get_cell, si, patch=(None,None,None)):
        if field.isUniform():
            cell_val = field.value()
            seq = patch2slice(patch)
            array[si, seq[0], seq[1], seq[2],:] = cell_val
        else:
            for cell_no in xrange(len(field)): 
                cell_coord = (float(cell_cx[cell_no]), float(cell_cy[cell_no]), float(cell_cz[cell_no]))
                cell = get_cell(cell_coord, dx, dy, dz)
                # Check if it is scalar
                if field_size == 1:
                    cell_val = field[cell_no]
                else:
                    cell_val = field[cell_no][:]
                # Add 1 due to low boundary
                array[si, cell[0], cell[1]+1, cell[2], :] = cell_val
    steps_index = {}
    for si, step in enumerate(time_steps):
        steps_index[step] = si
        time_dir = os.path.join(casedir, "%s/%s" % (timestep_foldername(step), field_name))
        field_file = ParsedParameterFile(time_dir)
        internal_field = field_file["internalField"]
        field2array(internal_field, get_cell_internal, si)
        if boundary:
            boundary_field = field_file["boundaryField"]
            for bname,face in boundary_patches:
                patch_field = boundary_field[bname]["value"]
                field2array(patch_field, get_cell_boundary, si, face)

    #NOTE: Only boundaries allowed in y for the moment
    x = np.linspace(bounds[0][0] + dx/2.0, bounds[1][0] - dx/2.0, ncx)
    z = np.linspace(bounds[0][2] + dz/2.0, bounds[1][2] - dz/2.0, ncz)
    if boundary:
        y = np.zeros(ncy_boundary)
        sly = slice(1, -1)
        y[0] = bounds[0][1]
        y[-1] = bounds[1][1]
        y[sly] = np.linspace(bounds[0][1] + dy/2.0, bounds[1][1] - dy/2.0, ncy)
    else:
        y = np.linspace(bounds[0][1] + dy/2.0, bounds[1][1] - dy/2.0, ncy)
        array = array[:,:,1:-1,:,:]

    return {"data": array, "axis": [x, y, z], "tidx": steps_index}
