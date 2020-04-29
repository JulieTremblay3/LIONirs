% VoxDims = [ Dimx, Dimy, Dimz ]
function VoxDims = get_PhysicalVoxDim( oMRI )
    VoxDims = [ oMRI.VoxDim.XColumn, oMRI.VoxDim.YSlice, oMRI.VoxDim.ZRow ];
 