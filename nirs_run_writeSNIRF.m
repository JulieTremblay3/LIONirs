function out = nirs_run_writeSNIRF(job)
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

%load(job.NIRSmat{1})     % Matrice NIRS
%Big loop over all subjects

    for Idx=1:numel(job.NIRSmat)   
        load(job.NIRSmat{Idx});
        nbch = numel(NIRS.Cf.H.C.n);
        dall = [];
        Sessionid = numel(NIRS.Dt.fir.pp); %prendre la dernière session
    try %create a defauld batchHistory file with lionirs pipeline used
    if isfield(job,'m_SNIRFBATCHhistory') %print session history a .m not a step could be save... 
        if job.m_SNIRFBATCHhistory
     if isfield(job,'c_SNIRFname')
            if isfield(job.c_SNIRFname,'b_SNIRFnamespecific')
                  name =  job.c_SNIRFname.b_SNIRFnamespecific.e_SNIRFnamespecific ;
                  outfilebatch= [job.f_writeNIRSdir{1},filesep,name,'BatchHistory.m']   ;     
            else
                outfilebatch= fullfile(job.f_writeNIRSdir{1},'BatchHistory.m');
            end
        else
            outfilebatch= fullfile(job.f_writeNIRSdir{1},'BatchHistory.m');
     end

    fid = fopen(outfilebatch,'w');
    for ibatch=1:numel(NIRS.Dt.fir.pp)
        jobinstruction = convertpretobatchname(NIRS.Dt.fir.pp(ibatch).pre);
        jobtext = ['matlabbatch{',num2str(ibatch),'}.',jobinstruction,'.'];       
        fieldlist = fieldnames(NIRS.Dt.fir.pp(ibatch).job);
        for idjob = 1:numel( fieldlist)
            id = 1;             
                if ~isstruct(eval(['NIRS.Dt.fir.pp(', num2str(ibatch),').job.' fieldlist{idjob}]))                                   
                    if iscell(eval(['NIRS.Dt.fir.pp(', num2str(ibatch),').job.' fieldlist{idjob}]))
                        logjob=    [jobtext, [fieldlist{idjob}],'={''',   eval(['NIRS.Dt.fir.pp(', num2str(ibatch),').job.' fieldlist{idjob}]),'''}'];
                        fprintf(fid,'%s',string(logjob));
                        fprintf(fid,'\n');
                        %sprintf('%s',string(logjob))
                    elseif isstr(eval(['NIRS.Dt.fir.pp(', num2str(ibatch),').job.' fieldlist{idjob}]))
                         logjob=    [jobtext, [ fieldlist{idjob}],'=''',   eval(['NIRS.Dt.fir.pp(', num2str(ibatch),').job.' fieldlist{idjob}]),''''];
                         fprintf(fid,'%s',string(logjob));   
                         fprintf(fid,'\n');
                         %sprintf('%sn',string(logjob))
                    elseif isnumeric(eval(['NIRS.Dt.fir.pp(', num2str(ibatch),').job.' fieldlist{idjob}]))
                        tmp = eval(['NIRS.Dt.fir.pp(', num2str(ibatch),').job.' fieldlist{idjob}]);
                        if numel(tmp) == 1 
                         logjob=    [jobtext, [fieldlist{idjob}],'=',   num2str(tmp)];
                         fprintf(fid,'%s',string(logjob));
                         fprintf(fid,'\n');
                        elseif size(tmp,1)<size(tmp,2)
                            logjob=    [jobtext, [fieldlist{idjob}],'=[',   num2str(tmp),']'];
                         fprintf(fid,'%s',string(logjob));
                         fprintf(fid,'\n');
                        elseif size(tmp,1)>size(tmp,2)
                             logjob=    [jobtext, [fieldlist{idjob}],'= [',   num2str(tmp'),']'];
                            fprintf(fid,'%s',string(logjob));
                            fprintf(fid,'\n');
                        end
                    end
                elseif isstruct(eval(['NIRS.Dt.fir.pp(', num2str(ibatch),').job.' fieldlist{idjob}]))
                    firstpart = ['NIRS.Dt.fir.pp(', num2str(ibatch),').job.' ];
                    temp =[fieldlist{idjob}];
                     fieldlist2 = fieldnames( eval([firstpart,temp]));
                     temp = [temp '.' fieldlist2{1}];
                    while isstruct(eval([firstpart,temp])) & isempty(isstruct(eval([firstpart,temp])));
                        fieldlist2 = fieldnames( eval([firstpart,temp]));
                        temp = [temp '.' fieldlist2{1}];                         
                    end
                           
                
                     if iscell(eval([firstpart,temp]))
                        logjob=    [jobtext, [temp],'={''',   eval([firstpart,temp]),'''}'];
                        fprintf(fid,'%s',string(logjob));
                        fprintf(fid,'\n');
                        %sprintf('%s',string(logjob))
                    elseif isstr(eval([firstpart,temp]))
                         logjob=    [jobtext, [temp],'=''',   eval([firstpart,temp]),''''];
                         fprintf(fid,'%s',string(logjob));   
                         fprintf(fid,'\n');
                         %sprintf('%sn',string(logjob))
                    elseif isnumeric(eval([firstpart,temp]))
                        tmp = eval(eval([firstpart,temp]));
                        if numel(tmp) == 1 
                         logjob=    [jobtext, [temp],'=',   num2str(tmp)];
                         fprintf(fid,'%s',string(logjob));
                         fprintf(fid,'\n');
                        elseif size(tmp,1)<size(tmp,2)
                            logjob=    [jobtext, [temp],'=[',   num2str(tmp),']'];
                         fprintf(fid,'%s',string(logjob));
                         fprintf(fid,'\n');
                        elseif size(tmp,1)>size(tmp,2)
                             logjob=    [jobtext, [temp],'= [',   num2str(tmp'),']'];
                            fprintf(fid,'%s',string(logjob));
                            fprintf(fid,'\n');
                        end                   
                     end
                    end
        end
    end           
     fclose(fid);
     disp(['Create log LIONirs as a Batch.m file : ', outfilebatch ])
        end
    end
        catch
            
        end

  
    for i = 1:numel(NIRS.Dt.fir.pp(1,Sessionid).p)
       [pathstr, name, ext] = fileparts(NIRS.Dt.fir.pp(1,Sessionid).p{i});
       d = fopen_NIR([pathstr,filesep,name,ext],nbch);       
        if ~isdir(job.f_writeNIRSdir{1})
           mkdir(job.f_writeNIRSdir{1});
       end
          
        if isfield(job,'c_SNIRFname')
            if isfield(job.c_SNIRFname,'b_SNIRFnamespecific')
                  name =  job.c_SNIRFname.b_SNIRFnamespecific.e_SNIRFnamespecific; 
                  outfile= [job.f_writeNIRSdir{1},filesep,name,'.snirf'];        
            else
                 if ~isdir([job.f_writeNIRSdir{1},filesep,name])
                     mkdir([job.f_writeNIRSdir{1},filesep,name]);
                 end 
                outfile= [job.f_writeNIRSdir{1},filesep,name,filesep,name,'.snirf'];
            end
        else
            outfile= [job.f_writeNIRSdir{1},filesep,name,filesep,name,'.snirf'];
             if ~isdir([job.f_writeNIRSdir{1},filesep,name])
                mkdir([job.f_writeNIRSdir{1},filesep,name]);
            end 
        end

       
        d = d'; %data time point x channels

            SD = [];
            SD.Lambda = NIRS.Cf.dev.wl;
            SD.SrcPos = NIRS.Cf.H.S.r.o.mm.p';
            SD.DetPos = NIRS.Cf.H.D.r.o.mm.p';
            

            SD.nSrcs = size(SD.SrcPos,1);
            SD.nDets = size(SD.DetPos,1);
            if numel(NIRS.Cf.H.C.id(2,:)) == size(d,2)
                ml = [NIRS.Cf.H.C.id(2,:)',NIRS.Cf.H.C.id(3,:)', ones(numel(NIRS.Cf.H.C.id(2,:)),1),NIRS.Cf.H.C.wl']; %srs,det,1,wav
                t = 1/NIRS.Cf.dev.fs:1/NIRS.Cf.dev.fs:size(d,1)*1/NIRS.Cf.dev.fs;
            end       
            
             aux5 =  NIRS.Dt.fir.aux5{1}; %trigger info to set in s structure
             %Example; [2, 9;
             %2, 1733;
             %1,1]
             %Trig 2 sample of time 2 et 1733, 
             %Trig 1 sample 
             nbmaxtrig =  max(aux5(:,1));
             s = zeros(numel(t), nbmaxtrig);
             for i = 1:size(aux5,1)
                 s(aux5(i,2),aux5(i,1))=1;
             end
                 
              SD.MeasList = ml;
            
             %need to have and aux file
            aux = zeros(size(d,1),1);
        if ~isfield(SD,'SpatialUnit')
            if mean(abs(SD.SrcPos(1,:)))>1 & mean(abs(SD.SrcPos(1,:)))<10 %probablement en cm
                SD.SpatialUnit = 'cm';  %or mm depending on your probe design
                disp(['Warning spatial unit field is missing , coordinate :',num2str(SD.SrcPos(1,:)),' seem to be in cm unit please verify your design.']);
            elseif mean(abs(SD.SrcPos(1,:)))>10 & mean(abs(SD.SrcPos(1,:)))<100 %probablement en mm
                SD.SpatialUnit = 'mm';  %or mm depending on your probe design
                disp(['Warning spatial unit field is missing , coordinate :',num2str(SD.SrcPos(1,:)),' seem to be in mm unit please verify your design.']);
            elseif mean(abs(SD.SrcPos(1,:)))<1
                SD.SpatialUnit = 'm';  %or mm depending on your probe design
                disp(['Warning spatial unit field is missing , coordinate :',num2str(SD.SrcPos(1,:)),' seem to be in meter unit please verify your design.']);
            end
        end
        try
            SD.Landmarks3D.pos = NIRS.Cf.H.F.r.o.mm;
            SD.Landmarks3D.labels =  {'Nz'        
                                    'LPA'
                                    'RPA'};
        catch %error compatibility mne si champs vide
            SD.Landmarks3D.pos = [1 0 0;
                                  0 -1 0;
                                  0 1 0];
            SD.Landmarks3D.labels =  {'Nz'        
                              'LPA'
                              'RPA'};
        end

        nirs= struct('d',d,'SD',SD,'t',t','s',s,'aux',aux);
        % Save SNIRF: Convert .nirs format data to SnirfClass object, save it to .snirf file (HDF5)
        fprintf('Saving %s ...\n', outfile);
        
   
        snirf_saved= SnirfClass(nirs);      
 
    %VERIFY IF 'ModifyBeerLambertLaw' WHERE APPLY
    %Change if unit DataType in measurementList when concentration when
    for i=1:numel(NIRS.Dt.fir.pp)
        if strcmp(NIRS.Dt.fir.pp(i).pre,'ModifyBeerLambertLaw')
          for k=1:numel(snirf_saved.data.measurementList) 
             if snirf_saved.data.measurementList(k).wavelengthIndex==1
                 SetDataType(snirf_saved.data.measurementList(k),1,'hbo'); 
             elseif snirf_saved.data.measurementList(k).wavelengthIndex==2
                 SetDataType(snirf_saved.data.measurementList(k),1,'hbr');
             end
          end
          disp('ModifyBeerLambertLaw applied data type is define as HbO and HbR');
        end
    end

        tic; snirf_saved.Save(outfile); toc  
        disp(['Save: ', outfile]);
    end
    end

    out.NIRSmat = job.NIRSmat;

end


function batchname = convertpretobatchname(pre)
  
    if strcmp(pre,'READ_RAW_NIRxScout') %READDATA
        batchname = 'spm.tools.nirsHSJ.M_readNIRS.E_readNIRxscout';
    elseif strcmp(pre,'READ_RAW_NIRSport'); 
        batchname =  'spm.tools.nirsHSJ.M_readNIRS.E_readNIRSport';
    elseif strcmp(pre,'READ_RAW_NIRS'); 
        batchname = 'spm.tools.nirsHSJ.M_readNIRS.E_rawhomer';
    elseif strcmp(pre,'READ SNIRF');
        batchname ='spm.tools.nirsHSJ.M_readNIRS.E_readSNIRF';
    elseif strcmp(pre,'readBOXY');
         batchname = 'spm.tools.nirsHSJ.M_readNIRS.boxy1';
    elseif strcmp(pre,'READ_RAW_BrainVision');
         batchname = 'spm.tools.nirsHSJ.M_readNIRS.E_GenericDataExportBV'; 
    elseif 0 %multimodal missing no pre... 
        batchname = 'spm.tools.nirsHSJ.M_readNIRS.M_readMultimodal.E_readEEG';          
    elseif strcmp(pre,'Concatenate nirs.mat Files')
        batchname = 'spm.tools.nirsHSJ.M_Segment.E_Concatenate_nirsmat'; %SEGMENT ONSET
    elseif strcmp(pre,'Concatenate File')
        batchname = 'spm.tools.nirsHSJ.M_Segment.E_Concatenate_file';   
    elseif strcmp(pre,'Segmentation')
        batchname =  'spm.tools.nirsHSJ.M_Segment.segment';
    elseif strcmp(pre,'Step Detection') %PREPROCESSING 
        batchname = 'spm.tools.nirsHSJ.M_preprocessing.E_artefactdetection';
    elseif strcmp(pre,'Filtered')
        batchname ='spm.tools.nirsHSJ.M_preprocessing.bpfilt';
    elseif strcmp(pre, 'Normalization')
        batchname = 'spm.tools.nirsHSJ.M_preprocessing.normalization';
    elseif strcmp(pre,'Nullify Bad Intervals')
        batchname = 'spm.tools.nirsHSJ.M_preprocessing.nullifybad';
    elseif strcmp(pre,'ModifyBeerLambertLaw')
        batchname ='spm.tools.nirsHSJ.M_preprocessing.ODtoHbOHbR';
    elseif strcmp(pre,'Epoch averaging')|strcmp(pre,'Epoch averaging (multiple files)')
         batchname ='spm.tools.nirsHSJ.M_preprocessing.E_average';
    elseif strcmp(pre, 'Prewhitened')    
         batchname ='spm.tools.nirsHSJ.M_preprocessing.E_prewhitening';
    elseif strcmp(pre, 'Detrend')  
         batchname ='spm.tools.nirsHSJ.M_preprocessing.E_detrend';
    elseif strcmp(pre, 'Cardiac')  %cardiac control pas de trace juste des canaux enlevé dans le nirs.mat 
         batchname ='spm.tools.nirsHSJ.M_preprocessing.E_chcardiaccontrol';
    elseif strcmp(pre,'Manual Gui')
        batchname =  'spm.tools.nirsHSJ.E_GUI';
    elseif strcmp(pre,'Manual Gui Segmentation')
        batchname =  'spm.tools.nirsHSJ.E_GUI';
    elseif strcmp(pre,'HyperScanCombine')
        batchname =  'pm.tools.nirsHSJ.M_Utility.M_HyperScanNIRS.E_HyperScanCombineNIRS';
    elseif strcmp(pre,'MarkCorrectedInYellow')
        batchname = 'spm.tools.nirsHSJ.M_dataComponent.E_MarkCorrectionAsNoise';

    

    elseif strcmp(pre,'Manual Gui Substract Component ')
         batchname ='spm.tools.nirsHSJ.M_dataComponent.E_substractcomponent';


    else
        batchname = pre;
    end
end