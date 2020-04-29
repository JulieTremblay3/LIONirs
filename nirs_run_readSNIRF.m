function out = nirs_run_readSNIRF(job)
%NIRS_RUN_READSNIRF This function takes in a job and creates an analysis
%step from it. It will read the snirf file and write a NIRS.mat construct
%that can be used in the next steps of an analysis.
%   prjname: .prj file for project
%   fileIn: raw data file in .snirf format
%   fileOut: .nir, .mat, .vhdr, .vmkr

% NOTE: While based on nirs_run_readNIRxscout, DATA.d in nirs_run_readSNIRF
% corresponds to the transposed of DATA.d in nirs_run_readNIRxscout. That
% means that while original code use the value of DATA.d as is, code reused
% from readNIRxscout has the DATA.d variable transposed (and the size
% function has the 1 argument switched to 2).

%PARAMETRE D'AQUISITION
NIRS.Cf.dev.n = 'SNIRF FILE';

    %Verifie si le path pointe vers un fichier
    pathout = job.output_path_SNIRF;
    if isfolder(pathout)
        choice = questdlg(['Directory ', pathout,' already exist do  you want to continue and overwrite this folder ?'], ...
        'WARNING',...
            'Yes Overwrite', 'Cancel','Cancel');
    switch choice
         case 'Cancel'
            return
        case 'Yes Overwrite'
            try
                 rmdir( pathout);
            catch
                msgbox(['Please remove ' , pathout,' manually'])
                return
            end
    end
    else
        mkdir(pathout)
    end

% MISC. DATA
% This section handle the data that is in the .snirf file but isn't part of
% the measurement per se. Such data include the subject age

inputsnirf = job.inputSNIRF;
assert(isfile(inputsnirf),"Input SNIRF doesn't exists."); % Check if the file still exists (could have been deleted in the meanwhile.)
[dir1,fil1,ext1] = fileparts(inputsnirf{1});
filehdr = fullfile(dir1,[fil1,ext1]);
importedData = H5ReadSNIRF(filehdr); % Imported data is the structure containing the information gathered from the SNIRF file by the H5ReadSNIRF function.

% Check if the coordinates are valid. Reject a coordinate if it is zero.
if ~(all(any(importedData.sourcePosition,2)) && all(any(importedData.detectorPosition,2)) && all(any(importedData.landmarkPosition,2)))
    error("Source, detector and landmark coordinates cannot be at the origin.");
end

% Check the spatial unit.
spatialUnitIndex = importedData.metaDataTags(:,1) == "SpatialUnit";
spatialUnit = importedData.metaDataTags(spatialUnitIndex,2);
switch spatialUnit
    case "mm"
        conversionRatio = 1; % 1 mm per mm.
    case "cm"
        conversionRatio = 10; % 10 cm per mm.
%   case "inch" % This may be a future feature if the .snirf specification adds more units.
    otherwise
        warning('No spatial units provided by SNIRF file.')
        conversionRatio = 1; % Assume the spatial unit to be mm.
end
% Enter source labels and position in the data struct.
NIRS.Cf.H.S.N = size(importedData.sourcePosition,1); %Nb of sources
NIRS.Cf.H.S.r.o.mm.p = importedData.sourcePosition' * conversionRatio; % Transpose the matrix to respect the convention and convert to mm.
NIRS.Cf.H.S.n = importedData.sourceLabels; % Import the source labels.

% Enter detector labels and positions in the data struct.
NIRS.Cf.H.D.N = size(importedData.detectorPosition,1); %Nb of detectors
NIRS.Cf.H.D.r.o.mm.p = importedData.detectorPosition' * conversionRatio; % Transpose the matrix to respect the convention and convert to mm.
NIRS.Cf.H.D.n = importedData.detectorLabels; % Import the detector labels.

%PARAMETRE
% Check if the subject age is supplied by the metadata.
ageIndex = importedData.metaDataTags(:,1) == "DateOfBirth";
subjectAge = importedData.metaDataTags(ageIndex,2);
if length(string(subjectAge)) == 1
    NIRS.Dt.s.age = subjectAge;  %job.age1 Add the age to the struct if the information is supplied.
else
    NIRS.Dt.s.age = ""; % Set to an empty string if the age is not available.
end
NIRS.Dt.ana.T1 = ''; %FICHIER DE IRM

%FILE
%NIRS.Dt.fir.pp.p %read and process snirf file

for Idx_File=1:numel(inputsnirf) % For each file in input.
    % DATA CHANNELS

%    if isfield(job.c_shortdistance,'b_shortdistanceavailable')
%        info.nbDet = size(ImportedData.detectorLabels,1);
%    end

    DATA.d = importedData.waveData; % DATA was part of the nirs_run_readNIRxscout function and was kept to avoid excessive refactoring. Could be replaced to improve performances and clarity.
    time = importedData.timeVector;
    if size(time,1) == 2 % If the time vector contains only two values, it implies that one is the initial time and the other the sampling interval.
        NIRS.Cf.dev.fs = 1/time(2); % Find the inverse of the sampling interval to get the sampling frequency.
    else
        assert(range(time(2:end)-time(1:end-1))==0,"The values are not evenly-spaced, the assumption concerning the sampling frequency does not hold.")
        NIRS.Cf.dev.fs = numel(time)/(time(end)-time(1)); % Otherwise, find the sampling frequency using the duration of the sampling and the number of data points. THIS ASSUMES EQUALLY-SPACED DATA POINTS.
    end
    NIRS.Cf.dev.wl = importedData.wavelengths;

    DATA.ml = [ [importedData.channel(:).sourceIndex]',[importedData.channel(:).detectorIndex]', zeros(length([importedData.channel(:).detectorIndex]),1), [importedData.channel(:).wavelengthIndex]' ]; % DATA.ml is an array containing all the information needed to describe a channel. Similar to the measurementList group of

    % Fiducial
    % NIRS.Cf.H.F.r.o.mm = sMtg.matFiducials(1:3,:)'; % Read fiducial information from the helmet montage file.
    NIRS.Cf.H.F.r.o.mm = importedData.landmarkPosition(:,1:3)'; % Read fiducial information from the SNIRF file.

% Handle response
    %ML
    %DISTANCE GEOMETRIQUE
    % Calculate probe distance (geometric distance between source and
    % detector) for each row of ml. Loop backward through the array so that
    % the memory space is preallocated.
    for channelIndex = size(DATA.d,1):-1:1
        temp{channelIndex,1} = ['S',num2str(DATA.ml(channelIndex,1)),'_D',num2str(DATA.ml(channelIndex,2))];  %Cell channel identification
        NIRS.Cf.H.C.id(1,channelIndex) = channelIndex;
        NIRS.Cf.H.C.id(2,channelIndex) = DATA.ml(channelIndex,1); %Source
        NIRS.Cf.H.C.id(3,channelIndex) = DATA.ml(channelIndex,2);  %et detecteur
        NIRS.Cf.H.C.wl(channelIndex) = DATA.ml(channelIndex,4);
        x1 =  NIRS.Cf.H.S.r.o.mm.p(1,DATA.ml(channelIndex,1));
        x2 =  NIRS.Cf.H.D.r.o.mm.p(1,DATA.ml(channelIndex,2));
        y1 =  NIRS.Cf.H.S.r.o.mm.p(2,DATA.ml(channelIndex,1));
        y2 =  NIRS.Cf.H.D.r.o.mm.p(2,DATA.ml(channelIndex,2));
        z1 =  NIRS.Cf.H.S.r.o.mm.p(3,DATA.ml(channelIndex,1));
        z2 =  NIRS.Cf.H.D.r.o.mm.p(3,DATA.ml(channelIndex,2));
        NIRS.Cf.H.C.gp(channelIndex,1) = sqrt((x2-x1)^2+(y2-y1)^2+ (z2-z1)^2);%À CALCULER DISTANCE GÉOMETRIQUE
        pos(channelIndex,1) = (x2-x1)/2+x1;
        pos(channelIndex,2) = (y2-y1)/2+y1;
        pos(channelIndex,3) = (z2-z1)/2+z1;
        pos(channelIndex,4) = NIRS.Cf.H.C.gp(channelIndex,1);
    end

    % Dist = NIRS.Cf.H.C.gp;
    %Enlever les canaux qui ne sont pas dans les distances
    %idbad = find(NIRS.Cf.H.C.gp > job.distmax  | NIRS.Cf.H.C.gp < job.distmin );
    %Enlever les canaux qui ne sont pas dans le mask

    % idbad = find([info.mask,info.mask ]==0);
    idbad = find(NIRS.Cf.H.C.gp<job.distmin|NIRS.Cf.H.C.gp>job.distmax);

    %ADD SHORT DISTANCE CHANNEL
    if isfield(job.c_shortdistance,'b_shortdistanceavailable')
        idbad = [idbad; find(DATA.ml(:,2)==16)];
        SSHORT = job.c_shortdistance.b_shortdistanceavailable.e_shortdistancesrs; %[2,3,5,7,11,10,13,15];
        tokeep = [];
        for isrs=1:numel(SSHORT)
         idkeep =  find(DATA.ml(idbad,1)==SSHORT(isrs)& DATA.ml(idbad,2)==16);
        % Keep all detector 8 source for short distance regression whenever the
        % distance on the helmet
         tokeep = [tokeep; idkeep];
        end
        idbad(tokeep) = [];
    end


     idall = ones(size(DATA.ml,1),1);
     idall(idbad) = [];
     %idall(idbad)=0;

     idgood = find(idall);
     for i = 1:numel(idgood)
        NIRS.Cf.H.C.n{i,1}  = temp{idgood(i) ,1};
     end
    NIRS.Cf.H.C.N = numel(NIRS.Cf.H.C.n);

   %  idbad = find(std(DATA.d)./mean(DATA.d)>job.m_removestd_tdms )
    NIRS.Cf.H.C.gp(idbad) = [];
    NIRS.Cf.H.C.wl(idbad) = [];
    NIRS.Cf.H.C.id(:,idbad)=[];
    DATA.ml(idbad,:)=[];
    try
    DATA.d(:,idbad)=[];
    catch
        msgbox('Do you record short distance probe ?')
    end
    %TRIG A VERIFIER ACTUELEMENT PRENDRE LE TRIG MANUELLE DE L'ORDI
    NIRS.Dt.fir.rons = [];

    NIRS.Cf.H.C.ok(:,Idx_File) = ones(size(DATA.ml,1),1);

   % NIRS.Cf.H.C.ok(find(stdbad),Idx_File)=0;
   % NIRS.Cf.H.C.ok(find(lowbad),Idx_File)=0;

    pos(idbad,:) = [];
    %Trigger value

    if isfield(job.c_shortdistance,'b_shortdistanceavailable')
        SHORTmlid = find(DATA.ml(:,2)== 16 & DATA.ml(:,4)== 1);
        distancemat = zeros(size(pos,1)/2,numel(SHORTmlid ));
        %DISTANCE BETWEEN CHANNEL FOR CLUSTERING
        for i= 1:size(pos,1)/2
            for j = 1:numel(SHORTmlid )
             x1 = pos(i,1);
             y1 = pos(i,2);
             z1 = pos(i,3);
             ishort = SHORTmlid (j);
             x2 = pos(ishort,1);
             y2 = pos(ishort,2);
             z2 = pos(ishort,3);
             distancemat(i,j) = sqrt((x2-x1)^2+(y2-y1)^2+ (z2-z1)^2);
            end
        end
    %take the closest from the regressor to create zone each channel is
    %take only once
       zone.SD.Lambda = NIRS.Cf.dev.wl;
       zone.SD.SrcPos = NIRS.Cf.H.S.r.o.mm.p';
       zone.SD.DetPos = NIRS.Cf.H.D.r.o.mm.p';
       zone.SD.nSrcs  = numel(zone.SD.SrcPos)/3;
       zone.SD.nDets  = numel(zone.SD.DetPos)/3;
       zone.SD.SrcAmp = ones(zone.SD.nSrcs ,2);
       zone.SD.DetAmp = ones(zone.SD.nDets,2);
       zone.color = jet(8);
       zone.ml = DATA.ml;
       zone.peak = [];
       zone.pos = pos;
       %identify the 8 first izone for each short channel position
       [tok,id] = sort(distancemat); %#ok<ASGLU>
       for izone=1:8
           id(izone, 1)
           zone.label{izone} = ['SHORT S', num2str(DATA.ml(id(1,izone),1))] ;
           zone.plot{izone} = [DATA.ml(id(1,izone),1),DATA.ml(id(1,izone),2)];
           zone.plotLst{izone} = id(1,izone);
       end
       %associate all other channel to the closest zone
       for ich = 1:size(distancemat,1)
          [val,izone]= min(distancemat(ich,:));
          if val>0
            zone.plot{izone} = [zone.plot{izone};[DATA.ml(ich,1),DATA.ml(ich,2)]];
            zone.plotLst{izone} = [zone.plotLst{izone}, ich];
          end

       end
        save( [dir1, '/SHORTDISTANCE.zone'],'zone');
    end

    % TRIGGERS / STIMULATION
    % The trigger/stimulation are events that are used to represent
    % specific stimuli (i.e. a task, or a sound) that may influence the
    % subject or perhaps only the start of the experiment.

    nbStim = length(importedData.stimulation); % Get the number of stimulations.

    evt = []; % Create an empty matrix for the trigger data.
    for k = nbStim:-1:1
        rawStim = importedData.stimulation(k).data; % Read the stimulation data from the .snirf file.
        if size(rawStim,2) == 1
           rawStim = rawStim'; % If the raw input is not wide enough, it means it has been encoded as a vector of dimension 3x1. As such, it need to be transposed.
        end
        assert(size(rawStim,2)==3,"Error. The size of the event data is not valid."); % Make sure that the value to append is correct.
        stimID = ones(size(rawStim,1),1); % Create a vector of the same height as the raw stim data. Overwrite the vector of the previous stimulation.
        trigIDnumber = str2double(importedData.stimulation(k).name); % Try to convert the stim name to a number.
        if isempty(trigIDnumber) % If it doesn't work, default to the stim group number. TODO Add a specific stim name to the NIRS.mat struct to conserve the name even if it is not numeric.
            trigIDnumber = k;
        end
        stimID(1:size(rawStim,1),1) = trigIDnumber; % Change every value of the vector for the name of the trig.
        evt = [evt ; [stimID,rawStim(:,1)]]; % Append the time/trigger pair for the current trigger to the event array.
    end
    evt = sortrows(evt,2); % Sort the array so that the triggers will be in chronological order.

    if numel(evt)>0 % Verify that the event matrix is not empty.
        NIRS.Dt.fir.aux5{Idx_File} = evt; % Assign the evt matrix to the data structure.
    else
        NIRS.Dt.fir.aux5{Idx_File} = [1 1]; % Create a row of two ones if no event is provided.
    end


        if ~isempty(evt)
            fid = fopen(fullfile(dir1,[fil1,'_evt1trig.m']),'w');
            fprintf(fid,'%s\n',['filename = [{''',fil1,'''}];' ]);
            fprintf(fid,'%s\n', 'Trigvalue = 1');
            if numel(evt)==1
                fprintf(fid,'%s\n', ['timingfile{1} = [',num2str(evt(i)/NIRS.Cf.dev.fs),'];']);
            else
                for i=1:numel(evt(:,1))
                    if i==1
                        fprintf(fid,'%s\n', ['timingfile{1} = [',num2str(evt(i,1)/NIRS.Cf.dev.fs),';']);
                    elseif i==numel(evt(:,1))
                        fprintf(fid,'%s\n', [num2str(evt(i,1)/NIRS.Cf.dev.fs),'];']);
                    else
                        fprintf(fid,'%s\n', [num2str(evt(i,1)/NIRS.Cf.dev.fs),';']);
                    end
                end
            end
            fclose(fid);
        end
         preprocess.STDamp =  job.STD_amp_choice.STD_amp;
         preprocess.STDampenable = job.STD_amp_choice.STD_enable;
         preprocess.STDmenu = job.STD_amp_choice.STD_menu;
         preprocess.DCamp = job.DC_amp_choice.DC_amp;
         preprocess.DCampenable = job.DC_amp_choice.DC_enable;
         fs=  NIRS.Cf.dev.fs;


        [dout,ok] = fast_preprocessing(DATA.d,  preprocess, fs); %#ok<ASGLU>
        NIRS.Cf.H.C.ok(:,Idx_File) = ok;
        NIRS.Dt.fir.sizebloc{Idx_File} = size(DATA.d,2); % 2 because it is transposed! Could also be size(DATA.d',1).
        %WRITE OUTPUT FILES
        NIRS.Dt.s.p = [pathout,'/'];%RAW NIRS PATH
        fileOut_nir = fullfile(pathout,[fil1,'.nir']);
        fileOutRoot_vhdr = fullfile(pathout,[fil1,'.vhdr']);
        fileOutRoot_vmrk = fullfile(pathout,[fil1,'.vmrk']);

        fwrite_NIR(fileOut_nir,DATA.d); % Untransposed.


        nirs_boxy_write_vhdr(fileOutRoot_vhdr,... %Output file
            fileOut_nir,... %DataFile
            fileOutRoot_vmrk,... %MarkerFile,...
            'nirs_convert_boxy',... %Function that created the header
            '',... %Channel Resolution
            '',... %Channel Units
            NIRS.Cf.H.C.n,... %names given as a column of cells
            1/NIRS.Cf.dev.fs*1e6,... %SamplingInterval in microseconds
            NIRS.Dt.fir.sizebloc{Idx_File}); %SamplingInterval in microseconds

        %dummy segment, add real marker read
        SD.Markers{1,Idx_File}.Type='New Segment';
        SD.Markers{1,Idx_File}.Description='';
        SD.Markers{1,Idx_File}.Position=1;
        SD.Markers{1,Idx_File}.Size=1;
        SD.Markers{1,Idx_File}.ChNumber=0; %all channels
        % temp_markers{1,1} = SD.Markers{1,Idx_File};
        idmarker = 1;
        for istimulus=1:size(evt,1)
            SD.Markers{idmarker,1}.Type='trigger';
            label = sprintf( 'S%3.0f', evt(istimulus, 1));
            SD.Markers{idmarker,1}.Description=label;
            SD.Markers{idmarker,1}.Position=evt(istimulus,2);
            SD.Markers{idmarker,1}.Size=1;
            SD.Markers{idmarker,1}.ChNumber=0; %all channels
            idmarker =  idmarker + 1;
        end
     %   temp_markers = SD.Markers{Idx_File,: };
        nirs_boxy_write_markers(fileOutRoot_vmrk,... %Output file
            fileOut_nir,... %DataFile
            SD.Markers);


        %FIRST STEP OF PROCESSING PP1
        NIRS.Dt.fir.pp(1).p{Idx_File,1} = fileOut_nir;
        NIRS.Dt.fir.pp(1).pre = 'READ_RAW_NIRxScout';
        NIRS.Dt.fir.pp(1).job = job;
end

        % Verify if the job creates or imports a project file.
        if isfield(job.c_createImportProjectSnirf, 'b_importProject')
            NIRS.Cf.H.prj = job.c_createImportProjectSnirf.b_importProject;
        else
            NIRS.Cf.H.prj = Create_PrjStruct(create_project_SNIRF(importedData)); % Create the project montage from the imported data.
        end

        % Save the project file.
        PrjData = NIRS.Cf.H.prj;
        [~, name, ~] = fileparts(inputsnirf{1});
        save([fullfile(pathout, name), '.prj'], '-struct', 'PrjData');

        % Save the NIRS.mat struct.
        save(fullfile(pathout,'NIRS.mat'),'NIRS');
        job.NIRSmat{1} = fullfile(pathout,'NIRS.mat');
        out.NIRSmat = job.NIRSmat;
end

function data = H5ReadSNIRF(filename)
    % This function reads a SNIRF file filename using high level H5 library
    % functions. It extracts the information of all the fields and returns
    % them as a structure.

    data = struct; % Initialise the return structure.

    % Read all nirs groups.
    nbNirs = H5CountGroups(filename, '/', '/nirs[0-9]*'); % Count the number of nirs groups (subjects) in the file.
    for nirsIndex = nbNirs:-1:1
        if nbNirs == 1 && H5CountGroups(filename, '/', '/nirs/') == 1
            nirsIndex = "";     %#ok<FXSET> % If there is only one nirs group, it is possible that there's no number added to it. Check and act accordingly.
        end                                 % This is probably a dirty hack, so watch out for unexpected behaviors.
        % Read channel info.
        dataInfo = h5info(filename, "/nirs" + nirsIndex + '/data1/dataTimeSeries');
        nbChannel = dataInfo.Dataspace.Size(1);
        data(nirsIndex).nbChannel = nbChannel;
        % Read informations of all channels by iterating through them.
        for n = nbChannel:-1:1
            infoChannelElement.sourceIndex = h5read(filename, "/nirs" + nirsIndex + "/data1/measurementList" + n + "/sourceIndex");
            infoChannelElement.detectorIndex = h5read(filename, "/nirs" + nirsIndex + "/data1/measurementList" + n + "/detectorIndex");
            infoChannelElement.wavelengthIndex = h5read(filename, "/nirs" + nirsIndex + "/data1/measurementList" + n + "/wavelengthIndex");
            infoChannelElement.dataType = h5read(filename, "/nirs" + nirsIndex + "/data1/measurementList" + n + "/dataType");
            infoChannelElement.dataTypeIndex = h5read(filename, "/nirs" + nirsIndex + "/data1/measurementList" + n + "/dataTypeIndex");
            data(nirsIndex).channel(n) = infoChannelElement;
        end

        % Read main data and time vector.
        data(nirsIndex).waveData = h5read(filename, "/nirs" + nirsIndex + "/data1/dataTimeSeries");
        data(nirsIndex).timeVector = h5read(filename, "/nirs" + nirsIndex + "/data1/time");

        % Read metadata.
        nbTags = H5CountGroups(filename, "/nirs" + nirsIndex, '/metaDataTags[0-9]*');
        for tagIndex = nbTags:-1:1
            data(nirsIndex).metaDataTags(tagIndex,1) = readUTF16Bytes(h5read(filename, "/nirs" + nirsIndex + "/metaDataTags" + tagIndex + "/name"));
            data(nirsIndex).metaDataTags(tagIndex,2) = readUTF16Bytes(h5read(filename, "/nirs" + nirsIndex + "/metaDataTags" + tagIndex + "/value"));
        end

        % Read probe information.

        % Surround label read code with try catch blocks as they are not
        % required by the format and may not be present. Assign them
        % default ("dummy") values if required.

        % Detector (Read position and labels).
        data(nirsIndex).detectorPosition = h5read(filename, "/nirs" + nirsIndex + '/probe/detectorPos');
        try
            for detectorIndex = size(data(nirsIndex).detectorPosition, 1):-1:1
                data(nirsIndex).detectorLabels(detectorIndex) = readUTF16Bytes(h5read(filename, "/nirs" + nirsIndex + '/probe/detectorLabels' + detectorIndex));
            end
        catch
            warning("No detector label detected, creating labels automatically.");
            data(nirsIndex).detectorLabels = strcat("D", string((1:size(data(nirsIndex).detectorPosition))')); % Create an array of numbered string of the length of detectorPosition. Ex: "D1"
        end
        % Source (Read position and labels).
        data(nirsIndex).sourcePosition = h5read(filename, "/nirs" + nirsIndex + '/probe/sourcePos');
        try
            for sourceIndex = size(data(nirsIndex).sourcePosition, 1):-1:1
                data(nirsIndex).sourceLabels(sourceIndex) = readUTF16Bytes(h5read(filename, "/nirs" + nirsIndex + '/probe/sourceLabels' + sourceIndex));
            end
        catch
            warning("No source label detected, creating labels automatically.");
            data(nirsIndex).sourceLabels = strcat("S", string((1:size(data(nirsIndex).sourcePosition))')); % Create an array of numbered string of the length of sourcePosition. Ex: "S1"
        end
        % Landmark (Read position and labels).
        data(nirsIndex).landmarkPosition = h5read(filename, "/nirs" + nirsIndex + '/probe/landmark');
        try
            for landmarkIndex = size(data(nirsIndex).landmarkPosition, 1):-1:1
                data(nirsIndex).landmarkLabels(landmarkIndex) = readUTF16Bytes(h5read(filename, "/nirs" + nirsIndex + '/probe/landmarkLabels' + landmarkIndex));
            end
        catch
            warning("No landmark label detected, creating labels automatically.");
            data(nirsIndex).landmarkLabels = strcat("L", string((1:size(data(nirsIndex).landmarkPosition))')); % Create an array of numbered string of the length of landmarkPosition. Ex: "L1"
        end
        % Wavelength
        data(nirsIndex).wavelengths = h5read(filename, "/nirs" + nirsIndex + '/probe/wavelengths');
        try
            data(nirsIndex).wavelengthsEmission = h5read(filename, "/nirs" + nirsIndex + '/probe/wavelengthsEmission');
        catch
            warning("No emission wavelength detected, defaulting to zero.");
            data(nirsIndex).wavelengthsEmission = 0;
        end

        % Read stimulation info.
        nbStim = H5CountGroups(filename, "/nirs" + nirsIndex + "/", '/stim[0-9]+'); % Count the number of stimulation groups in the file.
        % Iterate throught each stimulation group.
        for n = nbStim:-1:1
            infoStimElement.name = readUTF16Bytes(h5read(filename, "/nirs" + nirsIndex + "/stim" + n + "/name")); % Read the data from the file and put it in an array, itself in a structure.
            infoStimElement.data = h5read(filename, "/nirs" + nirsIndex + "/stim" + n + "/data");
            data(nirsIndex).stimulation(n) = infoStimElement; % Create an array and populate it.
        end
        if  isa(nirsIndex, 'string')
            break % Break from the loop if there is only unnumbered nirs group.
        end
    end
end

function numberOfGroups = H5CountGroups(filename, path, regEx)
    % This function counts the number of groups satisfying the regEx
    % pattern at the specified path of the specified HDF5 file.
    %   Ex: nb = H5CountGroups('C:\test\file.h', '/group1/subgroup2', '/word[0-9]+');
    pathInfo = h5info(filename, path);              % The the information of the specified path.
    groups = [pathInfo.Groups(:).Name];             % Get a character vector with the name of each group in path in the filename file.
    numberOfGroups = length(regexp(groups, regEx)); % Find the number of groups satisfying the pattern.
end

function str = readUTF16Bytes(bytesVector)
    % This function convert a vector of uint8 to characters and then
    % converts those characters to a string.
    % It uses US-ASCII, then UTF-16, depending on the number of word
    % characters detected after decoding.
    charVector = native2unicode(bytesVector, 'US-ASCII')';      % Try to read the byte array using US-ASCII encoding.
    if 0.90 > count(charVector,'\w')/length(charVector)         % If word characters make up less than 90% of the string...
        charVector = native2unicode(bytesVector, 'UTF-16')';    % Switch to UTF-16 encoding.
    end
    str = string(charVector);
end
