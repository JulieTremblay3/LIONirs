function [data,infoBV,label,ind_dur_ch] = fopen_EEG(varargin)
%fopen(filename) %open whole bloc
%fopen(filename tstart, tstop) open a segment only
%Read EEG.dat format analyser generic data export
%Input file EEG generic data export .dat, .vmrk, .vhdr
%Output data will be in vectorized (time x channel matrix) event if they
%are save in multiplex you could use the output of this function directly in vectorised
%data matrix for the whole file
%info (vhdr information)
%vmrk 
try 
  if numel(varargin) >= 1
    filename = varargin{1};
  end
  if numel(varargin) >= 2
    tstart = varargin{2};
  else
    tstart = 1;
  end
  if numel(varargin) >= 3
    tstop = varargin{3};
  end
  
    [pathdata, filenametot, ext]=fileparts(filename);
    [label,ind_dur_ch] = read_vmrk_all(fullfile(pathdata,[filenametot,'.vmrk' ]));
    %open a EEG data file and return the data
    infoBV = read_vhdr(fullfile(pathdata,[filenametot,'.vhdr' ]));
        fid=fopen(fullfile(pathdata,[filenametot,ext]));
        if isfield(infoBV,'DataFormat')
            if strcmp(infoBV.DataFormat,'ASCII')
                 data = [];
                if strcmp(infoBV.DataOrientation,'VECTORIZED')
                    while ~feof(fid)
                        tline = fgetl(fid);
                        [token,remain] =  strtok(tline,' ');
                        k = strfind(remain,',');
                        remain(k)='.';
                        data = [data,str2num(remain)'];
                    end
                elseif strcmp(infoBV.DataOrientation,'MULTIPLEXED')
                     tline = fgetl(fid); %first line to discart label of data
%                      while ~feof(fid)
%                          tline = fgetl(fid);
                        A = fscanf(fid,'%c');
                         k = strfind(A,','); %ensure comma tranfer in dot 
                         A(k)='.';
                         data = [data,str2num(A)];
                        
    %                 end
                end
                
                1
            elseif strcmp(infoBV.DataFormat,'BINARY')
                data= fread(fid,inf,'float32',0,'ieee-le');
                if strcmp(infoBV.DataOrientation,'VECTORIZED')
                    infoBV.DataPoints =numel(data)/infoBV.NumberOfChannels;
                    data = reshape(data, infoBV.DataPoints,infoBV.NumberOfChannels);
                elseif strcmp(infoBV.DataOrientation,'MULTIPLEXED')
                    infoBV.DataPoints =numel(data)/infoBV.NumberOfChannels;
                    data = reshape(data, infoBV.NumberOfChannels,infoBV.DataPoints)';
                end
            end
        end
    fclose(fid);
   
      if tstart<0
            padding = zeros(round(abs(tstart)./((infoBV.SamplingInterval)*1e-6)),size(data,2));
            data = [padding; data];
      end
    if numel(varargin)>2 
        pstart = round(tstart/((infoBV.SamplingInterval)*1e-6));
        if pstart <= 0
            pstart = 1;
        end
        pstop = round(tstop/((infoBV.SamplingInterval)*1e-6));
        tmax = size(data,1)*infoBV.SamplingInterval*1e-6; 
        if pstop > size(data,1)
           % data= [data;flipud(data)];
             data= [data;nan(size(data))];
            data = data(pstart:pstop,:);
            disp(['Warning out of range padding in the last bloc aux file: ', filename ' stop at ',num2str(tmax), ' and need to pad time until ',num2str(tstop)]);
        else
            data = data(pstart:pstop,:);
        end
        infoBV.DataPoints = size(data,1);
    end
        
catch
    disp(['Failed to open file ', filename,' use binary export format'] );
end
end
function infoBV = read_vhdr(file)
fid = fopen(file);

while  ~feof(fid)
    fline = fgetl(fid);
    %if ~ischar(fline),   break,   end
    if strcmp(fline,'[Channel Infos]')
        mode = 'Channelinfo';
        ilabel = 1;
    elseif strcmp(fline,'[Coordinates]')
        mode = 'Coordinates';
        ilabel = 1;
    elseif strcmp(fline,'[Comment]')
        break
    end    
    [rem,tok]=strtok(fline,'=');
    if numel(rem)<2
    elseif  strcmp(rem,'NumberOfChannels')
        infoBV.NumberOfChannels = str2num(tok(2:end));
    elseif strcmp(rem,'DataType')
        infoBV.DataType = tok(2:end);
    elseif strcmp(rem,'DataFormat')
        infoBV.DataFormat = (tok(2:end));
    elseif strcmp(rem,'DataOrientation')
        infoBV.DataOrientation = (tok(2:end)); 
    elseif strcmp(rem,'BinaryFormat')
        infoBV.BinaryFormat = tok(2:end);
    elseif strcmp(rem,'DataPoints')
        infoBV.DataPoints = str2num(tok(2:end));
    elseif strcmp(rem,'SamplingInterval')
        infoBV.SamplingInterval=str2num(tok(2:end));
    elseif strcmp(rem,'DecimalSymbol')
        infoBV.DecimalSymbol=str2num(tok(2:end)); 
    elseif strcmp(rem,'Layers')
        infoBV.Layers=str2num(tok(2:end));
    elseif strcmp(rem,'SegmentDataPoints')
        infoBV.SegmentDataPoints = str2num(tok(2:end));
    elseif strcmp(rem,'SegmentDataPointsPre')
        infoBV.SegmentDataPointsPre = str2num(tok(2:end));  
    elseif strcmp(rem,'SegmentDataPointsPost')
        infoBV.SegmentDataPointsPost = str2num(tok(2:end));  

    elseif strcmp(rem(1:2),'Ch') %& strcmp(tok(1:5),'=SEEG'
        switch mode
            case 'Channelinfo'
                [rem,tok]=strtok(tok(2:end),',');
                infoBV.name_ele{ilabel} = rem; %(7:12);
                ilabel = ilabel+1;
            case 'Coordinates'
                 [rem,tok]=strtok(tok(2:end),',');
                 infoBV.coor_r(ilabel) =str2num(rem);
                 [rem,tok]=strtok(tok,',');
                 infoBV.coor_theta(ilabel) = str2num(rem);
                 [rem,tok]=strtok(tok,',');
                 infoBV.coor_phi(ilabel) = str2num(rem);
                 ilabel = ilabel+1;
        end
    end
  
end
fclose(fid);
end