function out = nirs_run_writeSNIRF(job)
%NIRS_RUN_WRITESNIRF This function export the data from the current job to
% a Shared Near Infrared File Format Specification (.snirf)-compliant file.

% TODO Handle batch processing of multiple subject to a single file.
% Can be achieved by iterating using the dataGroupsNumber variable.

    % Verify what is the selected export mode.
    if isfield(job.B_SNIRFExportLocation,'B_SNIRFCreate') % Check if in Create mode.
        
        subjectIndex = 1; % Create a new file, so the subject is the first piece of data.
        
        % Verify if the folder already exists %
        pathout = job.B_SNIRFExportLocation.B_SNIRFCreate.SNIRF_outpath;
        if isfolder(pathout)
            choice = questdlg(['Directory ', pathout,' already exist do you want to continue and overwrite this folder ?'], ...
            'WARNING',...
                'Yes Overwrite', 'Cancel','Cancel');
            switch choice
                case 'Cancel'
                    return
                case 'Yes Overwrite'
                    try 
                        rmdir(pathout);
                    catch
                        msgbox(['Please remove ' , pathout,' manually'])
                    return
                    end
            end
        else
            mkdir(pathout)
        end

        % Filename of the created file.
        filename = fullfile(job.B_SNIRFExportLocation.B_SNIRFCreate.SNIRF_outpath, (job.B_SNIRFExportLocation.B_SNIRFCreate.FileOutput+".snirf") ); 
        
        % Check if the file already exist.
        if isfile(filename)
        choice = questdlg(['File ', filename,' already exist do you want to continue and overwrite this file ?'], ...
        'WARNING',...
            'Yes Overwrite', 'Cancel','Cancel');
            switch choice
                case 'Cancel'
                    return
                case 'Yes Overwrite'
                    try 
                        % delete(filename); % Would not be safe as the filename
                        % could contain wildcard characters.
                        warning('The file overwrite function is not enabled as the filename is not properly sanitized.');
                    catch
                        msgbox(['Please remove ' , filename,' manually'])
                        return
                    end
            end
            % Open the existing file
            fileID = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
        else
            % Create new file.
            fcpl = H5P.create('H5P_FILE_CREATE');
            fapl = H5P.create('H5P_FILE_ACCESS');
            fileID = H5F.create(filename,'H5F_ACC_TRUNC',fcpl,fapl);
        end
        
    elseif isfield(job.B_SNIRFExportLocation,'B_SNIRFAppend') % Check if in Append mode.
        % Check whether the selected file exists.
        filename = job.B_SNIRFExportLocation.B_SNIRFAppend.f_SNIRFfile;
        if isfile(filename)
            nbNirsGroups = H5CountGroups(filename,'/','/nirs[0-9]*');
            subjectIndex = nbNirsGroups + 1; % Set the nirsgroup number to the next available index. Should work as group number should be sequentials.
            fileID = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
        else
            error("The selected file does not exist.")
        end
    else
        error("Error. No valid option selected.")
    end

    % Prepare the required variables %
    NIRSTemp = load(job.NIRSmat{1}); % Import the NIRS matrix that contains the general information and links to the data.
    NIRS = NIRSTemp.NIRS;
    rDtpTemp = NIRS.Dt.fir.pp(1); % Data filepath.
    rDtp = rDtpTemp.p;
    NC = NIRS.Cf.H.C.N; % Channel name.
    waveData = fopen_NIR(rDtp{subjectIndex},NC); % Import data into the workspace.
    stimulationTemp = NIRS.Dt.fir.aux5;
    stimulation = stimulationTemp{subjectIndex};
    
    % Type of int used for indices in the .snirf file (must be the same).
    intType = 'H5T_NATIVE_INT'; % SNIRF accepts 8-bits, 16-bits, 32-bits and 64-bits integers.
    intConv = @(n) int32(n); % Function handle for the conversion of the index value to the chosen integer format.
    floatType = 'H5T_NATIVE_DOUBLE';
    
    % Datatype of the measurements.
    datatype = 001; % This corresponds to continuous wave amplitude measurement.
    datatypeIndex = 1; % Data-type specific parameter indices. A vector describing additional information depending on the datatype. Maps to data in the probe group.
    
    % Data used to create the file.
    formatVersion = "1.10";                     % Version of the .snirf format. The most recent version as of 22th of May 2019 is V1.10.
    dataGroupsNumber = 1;                       % Number of data groups. Used to iterate through the data groups of the SNIRF file.
    % sourceNumber = NIRS.Cf.H.S.N;             % Number of sources.
    % detectorNumber = NIRS.Cf.H.D.N;           % Number of detectors.
    % channelNumber = sourceNumber * detectorNumber;
    
    % Fiducial landmark data labels.
    landmarkLabels = ["NAS", "LPA", "RPA"]; % Hardcoded order in the mat structure.
    
    % MAIN FILE CREATION %
    
    % SNIRF files are HDF5 files. As such, they are constructed using
    % groups and datasets that act like folders and files.
    %       Ex: The index of the detector of the second channel of the
    %       first data group is located at :
    %           /data1/measurementList2/detectorIndex
    
    % Data field initialisation
    for groupIndex = dataGroupsNumber:-1:1
        
        % Create the nirs (subject group).
        h5createGroup(fileID, "/nirs" + subjectIndex);
        
        dataTimeSeries = waveData; % Put the data array in the dataTimeSeries array.
        
        samplingRate = NIRS.Cf.dev.fs; % Get the sampling rate.
        samplingInterval = 1/samplingRate; % Get the sampling interval.
        % time = (0:size(waveData,2))*samplingInterval; % Create a linear space of the same length as the data and with value at samplingInterval of eachothers.
        time = [0 ; samplingInterval]; % Second option for time data (preferable because the data is measured at a fixed frequency).
        
        % Create the data group.
        h5createGroup(fileID, "/nirs" + subjectIndex + "/data" + groupIndex);
        
        % Write the data array and the time vector.
        h5createwrite(fileID, "/nirs" + subjectIndex + "/data" + groupIndex + "/time", time, floatType); 
        h5createwrite(fileID, "/nirs" + subjectIndex + "/data" + groupIndex + "/dataTimeSeries", dataTimeSeries, floatType);
        
        % Channel information initialisation
        channelData = NIRS.Cf.H.C;
        for channelIndex = size(waveData,1):-1:1
            % Create the measurementList group.
            h5createGroup(fileID, "/nirs" + subjectIndex + "/data" + groupIndex + "/measurementList" + channelIndex);
            % Create the datasets.
            h5createwrite(fileID, "/nirs" + subjectIndex + "/data" + groupIndex + "/measurementList" + channelIndex + "/sourceIndex", intConv(channelData.id(2,channelIndex)), intType);
            h5createwrite(fileID, "/nirs" + subjectIndex + "/data" + groupIndex + "/measurementList" + channelIndex + "/detectorIndex", intConv(channelData.id(3,channelIndex)), intType);
            h5createwrite(fileID, "/nirs" + subjectIndex + "/data" + groupIndex + "/measurementList" + channelIndex + "/wavelengthIndex", intConv(channelData.wl(channelIndex)), intType);
            h5createwrite(fileID, "/nirs" + subjectIndex + "/data" + groupIndex + "/measurementList" + channelIndex + "/dataType", intConv(datatype), intType);
            h5createwrite(fileID, "/nirs" + subjectIndex + "/data" + groupIndex + "/measurementList" + channelIndex + "/dataTypeIndex", intConv(datatypeIndex), intType);
        end
    end
    
    % Initialise version
    h5CreateStringDataset(fileID, "/formatVersion", formatVersion);
    
    % Stim field initialisation
    uniqueStimulation = unique(stimulation(:,1)); % Find every unique stimulation identifier.
    
    for stimIndex = length(uniqueStimulation):-1:1                                  % Iterate through the unique stimulation identifiers.
        stimName = string(uniqueStimulation(stimIndex));                            % Set the name of the stimulation.
        stimArrayIndex = find(stimulation(:,1) == uniqueStimulation(stimIndex));    % Find the indices of every value that correspond to the stimulation identifier.
        stimElement = zeros(size(stimArrayIndex,1),3);                              % Create an array of the right size.
        stimElement(:,1) = stimulation(stimArrayIndex,2);                           % Find every value that correspond to the given indices and add them to the array.
        
        % stimDataElement(:,2) = 0; % Redundant as the array is created
        % filled with zeros already. This row represent the duration of the
        % stimuli.
        
        stimElement(:,3) = 1; % Fill the last collumn with 1 as the stimulation value is 1 (true/active).
        
        % Create the stimgroup.
        h5createGroup(fileID, "/nirs" + subjectIndex + "/stim" + stimIndex);
        % Write the stim data.
        h5createwrite(fileID, "/nirs" + subjectIndex + "/stim" + stimIndex + "/data", stimElement, floatType); % Write the stimulation data to the HDF5 file.
        h5CreateStringDataset(fileID, "/nirs" + subjectIndex + "/stim" + stimIndex + "/name", stimName); % Write the stimulation id to the file.
    end
    
    % Probe field initialisation
    wavelenghts = NIRS.Cf.dev.wl;
    sourcePos = NIRS.Cf.H.S.r.o.mm.p'; % The source and detector position have to be transposed to respect the specifications.
    detectorPos = NIRS.Cf.H.D.r.o.mm.p';
    % sourceLabels = NIRS.Cf.H.S.n; % Not in NIRS.mat despite nirsHSJ documentation.
    % detectorLabels = NIRS.Cf.H.D.n;
    landmark = NIRS.Cf.H.F.r.o.mm'; % Same as for the detector and source coordinates. Last collumn is filled with zeros because no label is defined.
    % landmarkLabels = NIRS.Cf.H.F.n; % Not in NIRS.mat
    useLocalIndex = 0;
    
    % Create the probe group.
    h5createGroup(fileID, "/nirs" + subjectIndex + "/probe");
    
    % Write probe information.
    h5createwrite(fileID, "/nirs" + subjectIndex + "/probe/wavelengths", wavelenghts, floatType);
    h5createwrite(fileID, "/nirs" + subjectIndex + "/probe/sourcePos", sourcePos, floatType);
    h5createwrite(fileID, "/nirs" + subjectIndex + "/probe/detectorPos", detectorPos, floatType);
    %{ 
    % Source labels (optional)
    for sourceLabelIndex = size(sourceLabels,1):-1:1
        h5CreateStringDataset(fileID, "/nirs" + subjectIndex + "/probe/sourceLabels" + sourceLabelIndex, sourceLabels(sourceLabelIndex));
    end
    %}
    %{
    % Detector labels (optional)
    for detectorLabelIndex = size(detectorLabels,1):-1:1
        h5CreateStringDataset(fileID, "/nirs" + subjectIndex + "/probe/detectorLabels" + detectorLabelIndex, detectorLabels(detectorLabelIndex));
    end
    %}
    h5createwrite(fileID, "/nirs" + subjectIndex + "/probe/landmark", landmark, floatType);
    % Landmark labels (optional)
    for landmarkLabelIndex = size(landmarkLabels,2):-1:1
        h5CreateStringDataset(fileID, "/nirs" + subjectIndex + "/probe/landmarkLabels" + landmarkLabelIndex, landmarkLabels(landmarkLabelIndex));
    end
    h5createwrite(fileID, "/nirs" + subjectIndex + "/probe/useLocalIndex", intConv(useLocalIndex), intType);
    
    % TODO Metadata initialisation
    metaDataTags = [ ["SubjectID", "TODO"]              % The "TODO" values
                     ["MeasurementDate", "TODO"]        % are placeholders
                     ["MeasurementTime", "TODO"]        % as the internal
                     ["SpatialUnit", "TODO" ]           % data format of
                     ["Model", string(NIRS.Cf.dev.n)]   % nirsHSJ doesn't
                     ["Age", string(NIRS.Dt.s.age)]     % yet support those
                     ["AcquisitionStartTime","TODO"]    % fields.
                     ["StudyID","TODO"]
                     ["StudyDescription","TODO"]
                     ["AccessionNumber","TODO"] ];
    for tagIndex = size(metaDataTags,1):-1:1
        h5createGroup(fileID, "/nirs" + subjectIndex + "/metaDataTags" + tagIndex); % Create the groupe and the datasets for each row of the array.
        h5CreateStringDataset(fileID, "/nirs" + subjectIndex + "/metaDataTags" + tagIndex + "/name", metaDataTags(tagIndex,1));
        h5CreateStringDataset(fileID, "/nirs" + subjectIndex + "/metaDataTags" + tagIndex + "/value", metaDataTags(tagIndex,1));
    end
    out = struct;
    out.NIRSmat = job.NIRSmat; % Pass the job to the next module.
end

    % Helper functions %
    
    function h5createGroup(fid,path)
        % This function creates a group named path in the file specified by
        % the fid.
        plist = 'H5P_DEFAULT';
        gid = H5G.create(fid,path,plist,plist,plist);
        H5G.close(gid);
    end
    
    function h5createwrite(fid,path,dataset,datatype)
        % This function creates a number array HDF5 dataset on the
        % desired path using the given dataset in file referenced by the
        % pid.
        
        % Find whether the array is a vector.
        if isvector(dataset) 
            % Set dimension.
            dimension = 1; % Must be made manually because for an array of dimension i.e. 1x5, ndims returns 2 while H5S.create_simple expects 1.
            % Set size, same as the input array.
            dataSize = length(dataset); % Cannot use size because it returns [1,n] but H5S.create_simple expects a scalar.
            % Set maximum dimension.
            maxdims = dataSize;
        elseif ismatrix(dataset)
            % Set dimension.
            dimension = 2;
            % Set size, same as the input array.
            dataSize = size(dataset);
            % Set maximum dimension.
            maxdims = size(dataset); % Create an array with the dimensions of the dataset.
        else
            error("Not a valid dimension.") % Raise an error if the dimensions cannot be handled by the function.
        end
        
        % IMPORTANT: The h5 library use row-major notation (C-like) for its
        % arrays while MatLab uses collumn-major (Fortran-like) notation:
        %   Ex:          [2,3]
        %       Row-major      Collumn-major
        %      +--+--+--+          +--+
        %      |  |  |  |          |  |
        %      +--+--+--+          +--+
        %                          |  |
        %                          +--+
        %                          |  |
        %                          +--+
        % That's why fliplr is required to convert the MatLab description
        % of the arrays' dimensions to a h5 compatible one.
        
        % Set h5 size.
        h5_dims = fliplr(dataSize);
        % Set maximum space.
        h5_maxdims = fliplr(maxdims);
        % Copy the propriety of the datatype.
        datatypeID = H5T.copy(datatype);
        % Create a dataspace for the array.
        dspace = H5S.create_simple(dimension,h5_dims,h5_maxdims); % Create a data space that has the dimension of the array.
        % Create dataset
        dset = H5D.create(fid,path,datatypeID,dspace,'H5P_DEFAULT'); % Create a new dataset.
        % Write data
        H5D.write(dset,datatypeID,'H5S_ALL','H5S_ALL','H5P_DEFAULT',dataset); % Write using the provided datatype and standard I/O.
        % Close file & resources
        H5T.close(datatypeID);
        H5S.close(dspace);
        H5D.close(dset);
    end

    function h5CreateStringDataset(fid,path,dataset)
        % This function converts a string scalar to a uint16 vector and
        % writes it to the HDF5 file.
        %       Important note, this function does not create a file, it
        %       needs to modify an already existing one.
        %       Derivative from code by: Jason Kaeding
        
        % Make sure the input is a string scalar.
        if ischar(dataset)
            dataset = string(dataset); % Convert character vectors to strings.
        end
        assert(isstring(dataset) && isscalar(dataset), "Error. The string input must be a string scalar.");
        
        % Convert the string to an uint16 array using UTF-16 encoding.
        bytesStr = unicode2native(dataset,'UTF-16');
        dataSize = length(bytesStr);
        maxSize = dataSize;
        
        % Set h5 size.
        h5_dims = dataSize;
        % Set maximum space.
        h5_maxdims = maxSize;
        % Copy the proprieties of the 8 bits datatype.
        datatypeID = H5T.copy('H5T_NATIVE_B8');
        % Create a dataspace for the array.
        dspace = H5S.create_simple(1,h5_dims,h5_maxdims); % Create a data space that has the dimension of the array.
        % Create dataset
        dset = H5D.create(fid,path,datatypeID,dspace,'H5P_DEFAULT'); % Create a new dataset.
        % Write data
        H5D.write(dset,datatypeID,'H5S_ALL','H5S_ALL','H5P_DEFAULT',bytesStr); % Write using the provided datatype and standard I/O.
        % Close file & resources
        H5T.close(datatypeID);
        H5S.close(dspace);
        H5D.close(dset);
    end

    function numberOfGroups = H5CountGroups(filename, path, regEx)
        % This function counts the number of groups satisfying the regEx
        % pattern at the specified path of the specified HDF5 file.
        %   Ex: nb = H5CountGroups('C:\test\file.h', '/group1/subgroup2', '/word[0-9]+');
        pathInfo = h5info(filename, path); % The the information of the specified path.
        groups = [pathInfo.Groups(:).Name]; % Get a character vector with the name of each group in path in the filename file.
        numberOfGroups = length(regexp(groups, regEx)); % Find the number of groups satisfying the pattern.
    end

