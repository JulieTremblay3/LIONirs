function Image = get_2DSegImage( oMRI, AxisPlane, PlaneNo )
    
    [nVoxX, nVoxY, nVoxZ] = size( oMRI.matSkinSegmentation );
    
    %  0----  X+   
    %  |   ________________          ________________   
    %  |  |                |        |                |   
    %  Z+ |                |        |                |   
    %     |                |        |                |   
    %     | Axial Slice    |        | Coronal Slice  |   
    %     | AxisPlane == 2 |     Y+ | AxisPlane == 3 |   
    %     | Through Y      |     |  | Through Z      |   
    %     |________________|     |  |________________|   
    %                            0---- X+              
    %      ________________
    %     |                | 
    %     |                |
    %     |                |
    %     | Sagittal Slice |
    %  Y+ | AxisPlane == 1 |   
    %  |  | Through X      |   
    %  |  |________________|
    %  0---- Z+        
    
    %Through X - Sagittal Slice
    if( AxisPlane == 1 )
        if( PlaneNo > nVoxX )
            Image = [];
            disp( 'Plane No Exceeds Image Dimension (ERR142)' );
            return;
        end
        Image = flipdim( reshape( oMRI.matSkinSegmentation( PlaneNo, :, : ), nVoxY, nVoxZ ), 1 );
    %Through Y - Axial Slice
    elseif( AxisPlane == 2 )
        if( PlaneNo > nVoxY )
            Image = [];
            disp( 'Plane No Exceeds Image Dimension (ERR142)' );
            return;
        end
        Image = reshape( oMRI.matSkinSegmentation( :, PlaneNo, : ), nVoxX, nVoxZ )';
    %Through Z - Coronal Slice
    elseif( AxisPlane == 3 )
        if( PlaneNo > nVoxZ )
            Image = [];
            disp( 'Plane No Exceeds Image Dimension (ERR142)' );
            return;
        end
        Image = flipdim( reshape( oMRI.matSkinSegmentation( :, :, PlaneNo ), nVoxX, nVoxY )', 1 );
    else
        disp( 'Incorrect Axis plane direction (ERR143)' );
    end