function newPrj=Createproject_forSDprobe(FileName)
%This function is to facilitate the project creation (.prj) file using file
%data already contain coordinate and source and detector position
%if the coordinate are in 2d visualisation will be over a flat 2d plane
%geometrie. If coordinate contain fiducial and are in 3d 
%FileName : .nirs with SD structure or snirf file 


[filepath,name,ext] = fileparts(FileName);
if strcmp(ext, '.nirs')
    NIRS= load(FileName,'-mat');
elseif strcmp(ext, '.snirf') 
    try
        snirf = SnirfClass(FileName);
    catch ME
        msgbox('Please install https://github.com/fNIRS/snirf_homer3 to support snirf class compatibility ')
        disp('Please install https://github.com/fNIRS/snirf_homer3 to support snirf class compatibility ')
        rethrow(ME) 
        return
    end
    NIRS.SD.Lambda = snirf.probe.wavelengths; 
    NIRS.SD.SrcPos = snirf.probe.sourcePos2D;
    NIRS.SD.DetPos = snirf.probe.detectorPos2D;
    NIRS.SD.SpatialUnit = 'mm';
end
%load NIRS get SD coordinate and montage information 
 newPrj = IO_Project_Data;
%load coordinate from NIRS file using SD.SrcPos SD.DetPos
handles.oDig_SubjectFiducials = Load_Digitization_SD(NIRS);
newPrj = set_Dig_SubjectFids( newPrj, handles.oDig_SubjectFiducials );
newPrj = set_Helmet( newPrj, FitHelmetOnSubject( get_Dig_SubjectFids(newPrj), get_Dig_CompleteHelmet(newPrj) ) );
newPrj = set_Helmet( newPrj, calc_Center2d( get_Helmet(newPrj) ) );

%set mtg here using ml structure
oHelmet = get_Helmet( newPrj);
new_Mtg = get_Mtg(oHelmet);
listsrs = [018001,020003,022005,024007,026009,028011,030013,032015,017002,019004,021006,023008,025010,027012,029014,031016,...
    050033,052035, 054037, 056039, 058041, 060043, 062045, 064047, 049034, 051036, 053038, 055040, 057042, 059044,061046, 063048,...
    082065,084067,86069,88071,90073,92075,94077,96079,81066,83068,85070,87072,89074,91076,93078,95080,...
    114097,116099,118101,120103,122105,124107,126109,128111,113098,115100,117102,119104,121106,123108,125110,127112 ];%Source POUR IDENTIFICATION USE IN HELMET HSJ PRJ STE-JUSTINE ISS
p  = 1;
new_Mtg.v_pSrc = [];
new_Mtg.v_pDet = [];
for i = 1:numel(NIRS.SD.SrcPos)/3
   new_Mtg.v_HolesMtg(p)= listsrs(i);
   new_Mtg.v_pSrc = [new_Mtg.v_pSrc,p];
   p = p+1; 
end
for i = 1:numel(NIRS.SD.DetPos)/3
   new_Mtg.v_HolesMtg(p)= i*1000000;
   new_Mtg.v_pDet = [new_Mtg.v_pDet,p];
   p = p+1;
end
new_Mtg.Gen_Params.Nb_Sources_Disp =  64;
new_Mtg.Gen_Params.Nb_Detect_Disp = 300;
new_Mtg.Gen_Params.AcqSystem = 'NIRx';
oHelmet = set_Mtg(oHelmet, new_Mtg);
newPrj = set_Helmet( newPrj, oHelmet );

%add flat2d surface as mri skin template file are located with the 3dMTG
%code
p = mfilename('fullpath');
[filepath,name,ext] =fileparts(p);

% VOXfile = fullfile([filepath,filesep, 'Template_Flat_2dsurface_fid'],'Flatd2fid.VOX')
% VOXhdrfile = fullfile([filepath,filesep, 'Template_Flat_2dsurface_fid'],'Flatd2fid.hdr')
% SRXfile = fullfile([filepath,filesep, 'Template_Flat_2dsurface_fid'],'Flat2d.SRX')
% SEGfile = fullfile([filepath,filesep, 'Template_Flat_2dsurface_fid'],'Flat2d.SEG')
% SEGhdrfile = fullfile([filepath,filesep, 'Template_Flat_2dsurface_fid'],'Flat2d.hdr')

%use MRI from previous build project 
% LoadedStruct = load(fullfile([filepath,filesep,'Template_Flat_2dsurface'],'2dFlatmontage.prj'),'-mat')
% prj_template = Load_PrjStruct( LoadedStruct, true );
% oMRI = get_MRI_Data(prj_template);
% newPrj = set_MRI_Data( newPrj, oMRI );