function display(obj)
%MRI_Data/DISPLAY command prompt call
disp(' ');
disp([inputname(1),' = ']);
disp( ' ' );
disp( 'MRI_Data:' );
disp( ' ' );
disp( sprintf(' size(matVoxels) = [%d,%d,%d]', size(obj.matVoxels,1), size(obj.matVoxels,2), size(obj.matVoxels,3) ) );
disp( sprintf('  dim(matVoxels) = [%d,%d,%d]', obj.VoxDim.XColumn, obj.VoxDim.YSlice , obj.VoxDim.ZRow ) );
disp( sprintf('     Voxel Depth = %d', obj.VoxDepth ) );
disp( sprintf('     Vox Max Int = %d', obj.Vox_IMax ) );
disp( sprintf('     Vox Min Int = %d', obj.Vox_IMin ) );
disp( sprintf('        Filename = %s', obj.PathFile ) );
disp( sprintf('      Nb Markers = %d', numel(obj.vMarkers) ) );
disp( sprintf('  Nb Vertex Skin = %d', numel(obj.Surfaces.Skin.VertexBuffer) ) );
disp( sprintf(' Nb Vertex CorHi = %d', numel(obj.Surfaces.CortexHiRes.VertexBuffer) ) );
disp( sprintf('Nb Vertex CorLow = %d', numel(obj.Surfaces.CortexLowRes.VertexBuffer) ) );
disp( ' Matrice de Fiducies =' );
disp( obj.matFiducials );
disp( ' Matrice de CoRegistration =' );
disp( obj.matCoRegistration );
disp( ' ' );