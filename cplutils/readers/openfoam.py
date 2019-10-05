import sys
import numpy as np
import glob
import os
import vtk
from vtk.util.numpy_support import vtk_to_numpy
# try:
#     from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile 
#     from PyFoam.RunDictionary.ParsedBlockMeshDict import ParsedBlockMeshDict
# except:
#     print "Error: PyFoam package is required"
#     sys.exit(1)

# def timestep_foldername(step):
#     ## NOTE: Only integer steps are supported for now
#     if step >= 1e6:
#         a = "{0:.16e}".format(step)
# 	print a.split('e')[0].rstrip('0').rstrip('.')+'e'+a.split('e')[1]
# 	return a.split('e')[0].rstrip('0').rstrip('.')+'e'+a.split('e')[1]
#     else:
#         return str(step)

def field_reader(casedir, field_name, time_steps=[], boundaries={}):
    def get_cell_data(reader, field):
        vtk_array = reader.GetOutput().GetCellData().GetArray(field)
        np_array = vtk_to_numpy(vtk_array)
        array_size = vtk_array.GetNumberOfComponents()
        return np_array, array_size

    def get_time_from_dirname(d):
        return float(os.path.basename(d).split("_")[1].split(".")[0])
    # Init regions present
    data = {"vol": {"data":[]}}
    for b in boundaries.keys():
        data[b] = {"data": []}
    for region in data.keys():
        if region is not "vol":
            patch = boundaries[region]
            casedir_vtk = os.path.join(casedir, "VTK/{}".format(patch))
        else:
            casedir_vtk = os.path.join(casedir, "VTK")

        if not os.path.exists(casedir_vtk):
            raise Exception("Directory '{}' not found. Use foamToVTK -useTimeName to generate the folder.".format(casedir_vtk))

        time_steps.sort() 
        time_steps_all = [get_time_from_dirname(d) for d in\
                        glob.glob(os.path.join(casedir_vtk, "*[0-9]*"))]
        folders_name = {get_time_from_dirname(d):d for d in\
                        glob.glob(os.path.join(casedir_vtk,"*[0-9]*"))}
        time_steps_all.sort()
        if not time_steps:
            time_steps = time_steps_all
        else:
            for s in time_steps:
                if s not in time_steps_all:
                    raise Exception("Timestep %s not present in openfoam case folder." % str(s))

        # Read time 0 for static data
        if region is not "vol":
            reader0 = vtk.vtkDataSetReader()
        else:
            reader0 = vtk.vtkUnstructuredGridReader()
        reader0.SetFileName(folders_name[0])
        reader0.Update()

        try:
            x, _ = get_cell_data(reader0, "ccx")
            x_axis = np.unique(x)
            ncx = len(x_axis)
            y, _ = get_cell_data(reader0, "ccy")
            y_axis = np.unique(y)
            ncy = len(y_axis)
            z, _ = get_cell_data(reader0, "ccz")
            z_axis = np.unique(z)
            ncz = len(z_axis)
            if region is not "vol":
                # NOTE: Assum patch perpendicular to -y direction
                ncells_patch = ncx*ncz
                cell_ids = range(0, ncells_patch)
            else:
                cell_ids, _ = get_cell_data(reader0, "cellID")
        except Exception:
            raise

        # Map coords to cell
        map_coord2cell = {}
        map_coord2patch_top = {}
        map_coord2patch_bot = {}
        for ix, xc in enumerate(x_axis):
            for iy, yc in enumerate(y_axis):
                for iz, zc in enumerate(z_axis):
                    map_coord2cell[(xc,yc,zc)] = (ix, iy, iz)

        field, field_size = get_cell_data(reader0, field_name)
        array = np.zeros((len(time_steps), ncx, ncy, ncz, field_size))
        steps_index = {}
        for si, step in enumerate(time_steps):
            steps_index[step] = si
            #
            reader = vtk.vtkDataSetReader()
            reader.SetFileName(folders_name[step])
            reader.Update()
            field, _ = get_cell_data(reader, field_name)
            # Internal cells
            for n, f in enumerate(field):
                cid = cell_ids[n]
                cell = map_coord2cell[(x[cid], y[cid], z[cid])]
                array[si, cell[0],cell[1],cell[2]] = f
        data[region]["data"] = array
        data[region]["axis"] = [x_axis, y_axis, z_axis]
        data[region]["tidx"] = steps_index
    # Concatenate arrays
    regions = ["vol"]
    if "bc" in data.keys():
        regions.insert(0,"bc") 
    if "wall" in data.keys():
        regions.append("wall")
    field_tuple = tuple([data[r]["data"] for r in regions])
    field_concat = np.concatenate(field_tuple, axis=2)
    yaxis_tuple = tuple([data[r]["axis"][1] for r in regions])
    yaxis_concat = np.concatenate(yaxis_tuple)
    tidx = data["vol"]["tidx"]
    return {"data": field_concat, "axis": [x_axis, yaxis_concat, z_axis], "tidx": steps_index}


#
# def field_reader(casedir, field_name, time_steps=[], boundary=False,
#                  boundary_patches=None):
#     if boundary:
#         if boundary_patches is None:
#             raise Exception("Boundary patches have to be provided.")
#     # Boundary fiel
#     time_steps.sort()
#     # time_steps_all = [int(float(os.path.basename(d))) for d in\
#     time_steps_all = [float(os.path.basename(d)) for d in\
#                       glob.glob(os.path.join(casedir,"[0-9]*"))]
#     folders_name = {float(os.path.basename(d)):d for d in\
#                       glob.glob(os.path.join(casedir,"[0-9]*"))}
#     print time_steps_all
#     time_steps_all.sort()
#     if not time_steps:
#         time_steps = time_steps_all
#     else:
#         for s in time_steps:
#             if s not in time_steps_all:
#                 print "Timestep %s not present in openfoam case folder." % str(s)
#                 raise Exception()
#     try:
#         cell_cx = ParsedParameterFile(os.path.join(casedir, "0/ccx"))["internalField"]
#         cell_cy = ParsedParameterFile(os.path.join(casedir, "0/ccy"))["internalField"]
#         cell_cz = ParsedParameterFile(os.path.join(casedir, "0/ccz"))["internalField"]
#     except Exception as e:
#         pass
#     try:
#         blockMesh = ParsedBlockMeshDict(os.path.join(casedir,
#                                         "constant/polyMesh/blockMeshDict"))
#     except Exception as e:
#         print e
#         pass
#
#     # Compute domain and mesh values
#     ncx, ncy, ncz = tuple(blockMesh["blocks"][2])
#     bounds = np.array(blockMesh.getBounds())
#     Lx, Ly, Lz = bounds[1] - bounds[0]
#     dx = Lx/ncx
#     dy = Ly/ncy
#     dz = Lz/ncz
#
#     # Get field dimension
#     time_dir = os.path.join(casedir, "%s/%s" % (str(0), field_name))
#     field = ParsedParameterFile(time_dir)["internalField"]
#     if type(field[0]) in [int, float]:
#         field_size = 1
#     else:
#         field_size = len(field[0])
#     #NOTE: Assumed in X and Z the boundaries are periodic
#     #NOTE: Top boundary = movingWall, Low boundary = CPLReceiveMD
#     #NOTE: Create an array assuming boundary=True, then slice at the end if False
#     ncy_boundary = ncy + 2
#     array = np.zeros((len(time_steps), ncx, ncy_boundary, ncz, field_size))
#
#     def get_cell_internal(cell_coord, dx,dy,dz):
#         return [int(cell_coord[0]/dx), int(cell_coord[1]/dy), int(cell_coord[2]/dz)]
#     def get_cell_boundary():
#         pass
#
#     def patch2slice(patch):
#         if patch == (None, None, None):
#             return (slice(None), slice(None), slice(None))
#         else:
#             return [int(-(c+1) / 2.0) if c != 0 else slice(None) for c in patch]
#
#
#     def field2array(field, get_cell, si, patch=(None,None,None)):
#         if field.isUniform():
#             cell_val = field.value()
#             seq = patch2slice(patch)
#             array[si, seq[0], seq[1], seq[2],:] = cell_val
#         else:
#             for cell_no in xrange(len(field)): 
#                 cell_coord = (float(cell_cx[cell_no]), float(cell_cy[cell_no]), float(cell_cz[cell_no]))
#                 cell = get_cell(cell_coord, dx, dy, dz)
#                 # Check if it is scalar
#                 if field_size == 1:
#                     cell_val = field[cell_no]
#                 else:
#                     cell_val = field[cell_no][:]
#                 # Add 1 due to low boundary
#                 array[si, cell[0], cell[1]+1, cell[2], :] = cell_val
#     steps_index = {}
#     for si, step in enumerate(time_steps):
#         steps_index[step] = si
#         time_dir = os.path.join(casedir, "%s/%s" % (folders_name[step], field_name))
#         field_file = ParsedParameterFile(time_dir)
#         internal_field = field_file["internalField"]
#         field2array(internal_field, get_cell_internal, si)
#         if boundary:
#             boundary_field = field_file["boundaryField"]
#             for bname,face in boundary_patches:
#                 patch_field = boundary_field[bname]["value"]
#                 field2array(patch_field, get_cell_boundary, si, face)
#
#     #NOTE: Only boundaries allowed in y for the moment
#     x = np.linspace(bounds[0][0] + dx/2.0, bounds[1][0] - dx/2.0, ncx)
#     z = np.linspace(bounds[0][2] + dz/2.0, bounds[1][2] - dz/2.0, ncz)
#     if boundary:
#         y = np.zeros(ncy_boundary)
#         sly = slice(1, -1)
#         y[0] = bounds[0][1]
#         y[-1] = bounds[1][1]
#         y[sly] = np.linspace(bounds[0][1] + dy/2.0, bounds[1][1] - dy/2.0, ncy)
#     else:
#         y = np.linspace(bounds[0][1] + dy/2.0, bounds[1][1] - dy/2.0, ncy)
#         array = array[:,:,1:-1,:,:]
#
#     return {"data": array, "axis": [x, y, z], "tidx": steps_index}
