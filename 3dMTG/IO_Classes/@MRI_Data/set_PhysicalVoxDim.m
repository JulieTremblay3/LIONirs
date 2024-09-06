% VoxDims = [ Dimx, Dimy, Dimz ]
function oMRI = set_PhysicalVoxDim( oMRI, VoxDims )
    oMRI.VoxDim.XColumn = VoxDims(1);
    oMRI.VoxDim.YSlice = VoxDims(2);
    oMRI.VoxDim.ZRow = VoxDims(3);
 