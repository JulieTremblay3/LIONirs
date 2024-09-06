path = 'G:\Fluence\MTG_MRI\'
%file = {'FVE27.prj',}
[file,path]=uigetfile('.prj','MultiSelect','on')
if ~iscell(file)
    file = {file}
end

%NORMALISATION 
%BATCH
for ifile=1:numel(file)
   Data_file =[path,file{ifile}];    
    LoadedStruct = load(Data_file,'-mat');
    PrjStruct = Load_PrjStruct(LoadedStruct,false);
    oHelmet = get_Helmet(PrjStruct);
    %Export elp and mtg
    Fileelp =[Data_file(1:end-4),'.elp'];
    Write_OPT3D_Elp( oHelmet, Fileelp,0,1);
    Filemtg = [Data_file(1:end-4),'.mtg'];
    save( [Filemtg ],'oHelmet');
end
% 
% newPrj = IO_Project_Data;
% handles.oDig_SubjectFiducials = Load_Digitization( Fileelp);
% newPrj = set_Dig_SubjectFids( newPrj, handles.oDig_SubjectFiducials );
% calc_Center( get_Helmet(newPrj));
% newPrj = set_Helmet( newPrj, calc_Center( get_Helmet(newPrj) ) );
% 
%  handles.oDig_MTGimport = load([path,Filemtg],'-mat');
%   if isfield(handles, 'oDig_MTGimport')
%         helmet1 = get_helmet(newPrj); 
%         helmet2 = handles.oDig_MTGimport.oHelmet;
%         helmet3 = combinehelmet(helmet1,helmet2);
%         newPrj = set_Helmet( newPrj,  helmet3 );
%   end
%   
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% FullFileName_Vox = 'G:\Fluence\MTG_MRI\NORS_NIRS\NORS.VOX'
% FullFileName_Hdr= 'G:\Fluence\MTG_MRI\NORS_NIRS\NORS.hdr'
%  fid_Hdr = fopen(FullFileName_Hdr);
%     if( fid_Hdr < 2 )
%         warndlg( 'Error loading Volume Data File (ERR113)' );
%     end
% 
%     %XRES
%     a_line=fgets(fid_Hdr);
%     if( a_line == -1 )
%         warndlg('Error loading Volume Data File (ERR114)');
%         return;
%     end
%     XRes=strread(a_line, '%d');
%     %set( handles.text_XRes, 'String', sprintf( '%d', XRes ) );
%     handles.XRes = XRes;
%                 
%     %YRES
%     a_line=fgets(fid_Hdr);
%     if( a_line == -1 )
%         warndlg('Error loading Volume Data File (ERR115)');
%         return;
%     end
%     YRes=strread(a_line, '%d');
%     %set( handles.text_YRes, 'String', sprintf( '%d', YRes ) );
%     handles.YRes = XRes;
%                 
%     %ZRES
%     a_line=fgets(fid_Hdr);
%     if( a_line == -1 )
%         warndlg('Error loading Volume Data File (ERR116)');
%         return;
%     end
%     ZRes=strread(a_line, '%d');
%     %set( handles.text_ZRes, 'String', sprintf( '%d', ZRes ) );
%     handles.ZRes = ZRes;
%                 
%     %Bytes per pixels (jumped)
%     a_line=fgets(fid_Hdr);
%     if( a_line == -1 )
%         warndlg('Error loading Volume Data File (ERR117)');
%         return;
%     end
%     BytesPerPix=strread(a_line, '%d');
%     handles.bpp = BytesPerPix;
%                 
%     %XZDim (vox dim inside a slice), in mm
%     a_line=fgets(fid_Hdr);
%     if( a_line == -1 )
%         warndlg('Error loading Volume Data File (ERR118)');
%         return;
%     end
%     XZVoxDim = strread(a_line);
%     if( numel(XZVoxDim) == 2 )
%         handles.XVoxDim = XZVoxDim(1);
%         handles.ZVoxDim = XZVoxDim(2);
%     else
%         handles.XVoxDim = XZVoxDim; %mm
%         handles.ZVoxDim = XZVoxDim; %mm
%     end
%     %set( handles.text_XDim, 'String', sprintf( '%1.4f', handles.XVoxDim ) );
%     %set( handles.text_ZDim, 'String', sprintf( '%1.4f', handles.ZVoxDim ) );
%     
%     %YDim (slice thickness) in mm
%     a_line=fgets(fid_Hdr);
%     if( a_line == -1 )
%         warndlg('Error loading Volume Data File (ERR119)');
%         return;
%     end
%     YDim=strread(a_line);
%     handles.YVoxDim = YDim; %mm
%     %set( handles.text_YDim, 'String', sprintf( '%1.4f', YDim ) );
% 
%     %Plane
%     a_line=fgets(fid_Hdr);
%     if( a_line == -1 )
%         warndlg('Error loading Volume Data File (ERR120)');
%         return;
%     end
%     Plane=strread(a_line, '%d');
%     %     if( Plane == 0 )
%     %         set( handles.text_Plane, 'String', 'Axial' );
%     %     elseif( Plane == 1 )
%     %         set( handles.text_Plane, 'String', 'Sagittal' );
%     %     elseif( Plane == 2 )
%     %         set( handles.text_Plane, 'String', 'Coronal' );
%     %     end
%     handles.Plane = Plane;
% 
%     %Jump this line - Skin treshold
%     a_line=fgets(fid_Hdr);
%     if( a_line == -1 )
%         warndlg('Error loading Volume Data File (ERR121)');
%         return;
%     end
% 
%     %Jump this line - Cortex treshold
%     a_line=fgets(fid_Hdr);
%     if( a_line == -1 )
%         warndlg('Error loading Volume Data File (ERR122)');
%         return;
%     end
%     
%     %LPA
%     a_line=fgets(fid_Hdr);
%     if( a_line == -1 )
%         warndlg('Error loading Volume Data File (ERR123)');
%         return;
%     end
%     Fids(2,:)=strread(a_line);
%     %if( numel( Fids(1,:) ) == 3 )
%     %    set( handles.text_LPA, 'String', sprintf('[%4d %4d %4d]', Fids(2,1), Fids(2,2), Fids(2,3) ) );
%     %end
%                 
%     %RPA
%     a_line=fgets(fid_Hdr);
%     if( a_line == -1 )
%         warndlg('Error loading Volume Data File (ERR124)');
%         return;
%     end
%     Fids(3,:)=strread(a_line);
%     %if( numel( Fids(2,:) ) == 3 )
%     %    set( handles.text_RPA, 'String', sprintf('[%4d %4d %4d]', Fids(3,1), Fids(3,2), Fids(3,3) ) );
%     %end
% 	
%     %NAS
%     a_line=fgets(fid_Hdr);
%     if( a_line == -1 )
%         warndlg('Error loading Volume Data File (ERR125)');
%         return;
%     end
%     Fids(1,:)=strread(a_line);
%     %if( numel( Fids(3,:) ) == 3 )
%     %    set( handles.text_NAS, 'String', sprintf('[%4d %4d %4d]', Fids(1,1), Fids(1,2), Fids(1,3) ) );
%     %end                
% 	
%     %Jump this line (VERTEX)
%     a_line=fgets(fid_Hdr);
%     if( a_line == -1 )
%         warndlg('Error loading Volume Data File (ERR126)');
%         return;
%     end
% 	
%     %Jump this line (ANTCOMM? - TAILARACH SYSTEM ORIGIN)
%     a_line=fgets(fid_Hdr);
%     if( a_line == -1 )
%         warndlg('Error loading Volume Data File (ERR127)');
%         return;
%     end
% 	
%     %Jump this line - ATLAS filename associated
%     a_line=fgets(fid_Hdr);
%     if( a_line == -1 )
%         warndlg('Error loading Volume Data File (ERR128)');
%         return;
%     end
% 	
%     %Jump this line - MIP filename associated
%     a_line=fgets(fid_Hdr);
%     if( a_line == -1 )
%         warndlg('Error loading Volume Data File (ERR129)');
%         return;
%     end
% 	
%     %Jump this line (INI)
%     a_line=fgets(fid_Hdr);
%     if( a_line == -1 )
%         warndlg('Error loading Volume Data File (ERR130)');
%         return;
%     end
%                 
%     %NOS
%     a_line=fgets(fid_Hdr);
%     if( a_line == -1 )
%         warndlg('Error loading Volume Data File (ERR131)');
%         return;
%     end
%     Fids(4,:)=strread(a_line);
%     %if( numel( Fids(4,:) ) == 3 )
%     %    set( handles.text_NOS, 'String', sprintf('[%4d %4d %4d]', Fids(4,1), Fids(4,2), Fids(4,3) ) );
%     %end
%     if isempty(find(Fids))     
%         warndlg('No fiducial point marked (ERR132)');
%         return;
%     end
%     handles.VOXFids = Fids;
%     
%     %----------------------------------------------------------------------
%     %--------------- SECTION VOX (BINARY) ---------------------------------
%     %----------------------------------------------------------------------
%     fid_Vox=fopen(FullFileName_Vox);
%     
%     if( BytesPerPix == 1 )
%         sDataType = 'uint8';
%     elseif( BytesPerPix == 2 )
%         sDataType = 'uint16';
%     else
%         warndlg('Error loading Volume Data File: unsupported bbp (ERR132)');
%         return;
%     end
% 
%     %Initialisation
%     countTotal = 0;
%     MaxPacketSize = 2500000;
%     VoxelsToRead = XRes*YRes*ZRes;
%     Data1D = zeros(VoxelsToRead, 1, sDataType);
%         
%     hWaitbar = waitbar(0,'Loading data...');
%     while( countTotal < VoxelsToRead )
%         
%         %Lit en format 'double': convertir en uint8/16
%         [Packet, count] =  fread(fid_Vox, MaxPacketSize, sDataType);
%         
%         if( BytesPerPix == 1 )% 8Bits
%             Data1D(countTotal+1:countTotal+count) = uint8(Packet);
%         else% 16Bits
%             Data1D(countTotal+1:countTotal+count) = uint16(Packet);
%         end
%         
%         countTotal = countTotal+count;
%         waitbar(countTotal/VoxelsToRead);
%     end
%     close(hWaitbar);
%     
%     if( BytesPerPix == 1 )
%         handles.matVoxels = uint8(reshape( Data1D, XRes, ZRes, YRes ));
%     else
%         handles.matVoxels = uint16(reshape( Data1D, XRes, ZRes, YRes ));
%     end    
%     fclose(fid_Hdr);
%     fclose(fid_Vox);
%     
% %     %%%%%%%%%%%%%%%%%%%%%%%%%
% % 'G:\Fluence\MTG_MRI\NORS_NIRS\cortex.SRX'
% % 'G:\Fluence\MTG_MRI\NORS_NIRS\cortex.SRX'
% % 'G:\Fluence\MTG_MRI\NORS_NIRS\skinlow.SRX'
% % 'G:\Fluence\MTG_MRI\NORS_NIRS\cortex.SEG'
% % 'G:\Fluence\MTG_MRI\NORS_NIRS\skin.SEG'
% % 'G:\Fluence\MTG_MRI\NORS_NIRS\brodmann.img'
% 
%   