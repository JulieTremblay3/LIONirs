function [stat neiglinkmat]= FindClusterBasedPermutationInMatrix(chanpos,neighbourdist,clustercritval, statobs, statrand, minnbchan  )
% Description: FindClusterBasedPermutationInMatrix
% define the neighboring rules to form a cluster in the case of a matrix 
% use clusterstat function to correct from multiple comparaison in the
% fieldtrip
% chanpos: channel position to define neighbour
% neighbourdist: Distance to be consider a neighbour
% clustercritval: Threshold to be significant
% statobs: stat observe
% statrand: stat random



if isempty(chanpos)
%     load('C:\data\EEGClinique\Script_Janie\chanpos.mat')
%    % neighbourdist='150';
%     disp('LOAD POUR JANIE')
end
for i=1:size(chanpos,1)
   listelectrode{i,1} = ['C',num2str(i)];    
end 
%neighbourdist = 'sig';
if strcmp(neighbourdist,'all')
   neiglinkmat(:,:)=ones(numel(statobs));  
   disp(['Neighborgs use all link'])
elseif strcmp(neighbourdist,'link')
    id = 1;
    nele = numel(listelectrode);
    matid = zeros(nele,nele);
    for ielex=1:nele
        ielex;
        ieley = 1;
        while ieley < ielex 
            option{id}.ele1 =  listelectrode{ielex};
            option{id}.ele2 =listelectrode{ieley};
            option{id}.matposition = [ielex,ieley];      
            matid(ielex,ieley)=id;
            neiglink(id).label1 = [listelectrode{ielex}];
            neiglink(id).label2 =  [listelectrode{ieley}];
             iele1(id)=ielex;
             iele2(id) =ieley;
%             neiglink(id).neighblabel =[neighbours(ielex).neighblabel;neighbours(ieley).neighblabel];
%             [seln] = match_str(listelectrode, neiglink(id).neighblabel);
%             neiglink(id).neigiele = seln;
            ieley = ieley + 1;
            id = id + 1;
        end
    end
    idhalf=find(matid);
    loc=matid(idhalf);
    %check link relate electrode
    neiglinkmat = false(length(neiglink),length(neiglink));

    %reoder with link
    iele_row_link =   iele1(loc);
    iele_col_link =   iele2(loc);

    for ilink=1:length(neiglink)
         listlink= find(iele_row_link  == iele_row_link(ilink) |iele_row_link   == iele_col_link(ilink)|iele_col_link   == iele_row_link(ilink)|iele_col_link   == iele_col_link(ilink) );
         neiglinkmat(ilink,listlink) = true;
    end   
   %  figure;imagesc(neiglinkmat)
       disp(['Neighborgs with commun link between them'])
elseif strcmp(neighbourdist,'sig') %allow neigborg among sig value
    neiglinkmat(:,:)=zeros(numel(statobs)); 
    idpos = find(statobs>clustercritval);
   % idminus = find(statobs<(-clustercritval));   
    for id = 1:numel(idpos)
        neiglinkmat(idpos(id),idpos)=1;
    end
%     for id = 1:numel( idminus)
%         neiglinkmat( idminus(id), idminus)=1;
%     end
    
elseif 0 %gere tout les liens du voisinage 

    neighbourdistnb = 4
    neighbours = compneighbstructfromgradelec(chanpos,listelectrode,neighbourdistnb);
   % neighbours = compneighbstructfromgradelec(chanpos,listelectrode,0);
    id = 1;
    nele = numel(listelectrode);
    matid = zeros(nele,nele);
    for ielex=1:nele
        ielex;
        ieley = 1;
        while ieley < ielex 
            option{id}.ele1 =  listelectrode{ielex};
            option{id}.ele2 =listelectrode{ieley};
            option{id}.matposition = [ielex,ieley];      
            matid(ielex,ieley)=id;
            neiglink(id).label1 = [listelectrode{ielex}];
            neiglink(id).label2 =  [listelectrode{ieley}];
             iele1(id)=ielex;
             iele2(id) =ieley;
            neiglink(id).neighblabel =[neighbours(ielex).neighblabel;neighbours(ieley).neighblabel];
            [seln] = match_str(listelectrode, neiglink(id).neighblabel);
            neiglink(id).neigiele = seln;
            ieley = ieley + 1;
            id = id + 1;
        end
    end
    idhalf=find(matid);
    loc=matid(idhalf);
    %check link relate electrode
    neiglinkmat = false(length(neiglink),length(neiglink));
    allneighb = {neiglink.neighblabel};
    %reoder with link
    iele_row_link =   iele1(loc);
    iele_col_link =   iele2(loc);
    %iele1 valeur de l'électrode 1 pour se lien (1 to 99)
    %iele2 valeur de l'électrode 2 pour se lien (1 to 99)
    %eiglink(ilink).neigiele liste des voisins associés au 2 electrode du lien.
    for ilink=1:length(neiglink)
        idloc = loc(ilink);
        %look for close neighbours one or the other of the electrode imply 
        listlink= find(sum([iele_row_link  == neiglink(idloc).neigiele ; iele_col_link  == neiglink(idloc).neigiele])>1);
        neiglinkmat(ilink,listlink) = true;       
    end 
    figure; imagesc( neiglinkmat)
 %   figure;imagesc(neiglinkmat)
    disp(['Neighborgs for clustering use optode with a distance ',num2str(neighbourdist),' between them'])
        
else %gere les combinaison commune 
      %distance entre les électrodes proche si [C1,C2,C3,C5]  et [C11,C12,c14,C15] sont voisin spatialement sont les liens considéré comme voisin entre les deux régions inclura 
    % les combinaisons des deux régions composés seront considéré.
    % C1 C11; C1 C12; C1 C14;  C1 C15;
    % C2 C11; C2 C12; C2 C14;  C2 C15; 
    % C3 C11; C3 C12; C3 C14;  C3 C15; 
    % C5 C11; C5 C12; C5 C14;  C5 C15;
  
   neighbourdistnb = str2num(neighbourdist);
   neighbours = compneighbstructfromgradelec(chanpos,listelectrode,neighbourdistnb);
     id = 1;
    nele = numel(listelectrode);
    matid = zeros(nele,nele);
   
         id = 1;
        for ieley=1:nele                 
              for  ielex=(ieley+1):nele
                 option{id}.ele1 =  listelectrode{ielex};
                option{id}.ele2 =  listelectrode{ieley};
                 option{id}.matposition = [ielex,ieley];      
                matid(ielex,ieley)=id;
                neiglink(id).label1 = [listelectrode{ielex}];
                neiglink(id).label2 =  [listelectrode{ieley}];
                idele1(id) = ielex;
                idele2(id) = ieley;
                label1{id} = listelectrode{ielex};
                label2{id} = listelectrode{ieley};   
                 neiglink(id).neighblabel1 = [neighbours(ielex).neighblabel];
                neiglink(id).neighblabel2 =   [neighbours(ieley).neighblabel];
                neiglink(id).ideleneigh1 =[neighbours(ielex).idneigh];
                neiglink(id).ideleneigh2=[neighbours(ieley).idneigh];
               id = id + 1;
             end
        
    end
    
    
    for  idhalf=find(matid);
        loc=matid(idhalf); 
        %check link relate electrode
        neiglinkmat = false(length(neiglink),length(neiglink));
    end
    
    if 0
     for ilink=1:length(neiglink)

        %look for close neighbours combinaisaion 
         %imply link 1
         tic
        [sellink1] = [match_str(label1,  label1(ilink)); match_str(label2,  label1(ilink))];
        [sellink2] = [match_str(label1,  label2(ilink)); match_str(label2,  label2(ilink))];
        temp = zeros(size(neiglinkmat,1),2);
        temp(sellink1,1)=1;
        temp(sellink2,2)=1;
         idfind=  find(sum(temp,2)==2);
         toc
        % disp(['find' label1(ilink), label2(ilink)])
           %idself
           for ineigh2 = 1:numel(neiglink(ilink).neighblabel2)
              [sellink1] = [match_str(label1,  label1(ilink)); match_str(label2,  label1(ilink))]; 
               label2neig = neiglink(ilink).neighblabel2{ineigh2 };
              [sellink2] = [match_str(label1,  label2neig); match_str(label2,  label2neig)];
               temp = zeros(size(neiglinkmat,1),2);
               temp(sellink1,1)=1;
               temp(sellink2,2)=1;
               idfind=  find(sum(temp,2)==2);
             disp(['find' label1(find(sum(temp,2)==2)), label2(find(sum(temp,2)==2))]);
               neiglinkmat(find(sum(temp,2)==2),find(sum(temp,2)==2)) = true;       
           end
            for ineigh1 = 1:numel(neiglink(ilink).neighblabel1)
              [sellink1] = [match_str(label1,  label2(ilink)); match_str(label2,  label2(ilink))]; 
               label1neig = neiglink(ilink).neighblabel1{ineigh1 };
              [sellink2] = [match_str(label1,  label1neig); match_str(label2,  label1neig)];
               temp = zeros(size(neiglinkmat,1),2);
               temp(sellink1,1)=1;
               temp(sellink2,2)=1;
               idfind=  find(sum(temp,2)==2);
               disp(['find' label1(find(sum(temp,2)==2)), label2(find(sum(temp,2)==2))]);
               neiglinkmat(find(sum(temp,2)==2),find(sum(temp,2)==2)) = true;       
           end
     end 
 
    
    end
     tic
    for  ilink=1:length(neiglink)
        %commun link
        idfind = find([(idele1(ilink)==idele1)|(idele1(ilink)==idele2) ] & [(idele2(ilink)==idele1)|(idele2(ilink)==idele2) ] );
       if ilink==3000
       1
       end
         neiglinkmat( ilink,idfind)= true;
          for ineigh2 = 1:numel(neiglink(ilink).neighblabel2)
              if idele1(ilink)~=neiglink(ilink).ideleneigh2(ineigh2)
               idfind = find([(idele1(ilink)==idele1)|(idele1(ilink)==idele2) ] & [(neiglink(ilink).ideleneigh2(ineigh2)==idele1)|(neiglink(ilink).ideleneigh2(ineigh2)==idele2) ]); 
               neiglinkmat( ilink,idfind)= true;
               %disp([idele1(idfind), idele2( idfind)])
               %disp([neiglink(idfind).label1,neiglink(idfind).label2 ]);
              end
          end
          for   ineigh1 = 1:numel(neiglink(ilink).neighblabel1)
              if  idele2(ilink)~=neiglink(ilink).ideleneigh1(ineigh1)
              idfind = find([(idele2(ilink)==idele1)|(idele2(ilink)==idele2) ] & [(neiglink(ilink).ideleneigh1(ineigh1)==idele1)|(neiglink(ilink).ideleneigh1(ineigh1)==idele2) ]); 
              neiglinkmat( ilink,idfind)= true;
               %disp([neiglink(idfind).label1,neiglink(idfind).label2 ]);
              %   disp([idele1(idfind), idele2( idfind)])
              end
          end
    end
    toc
end  
   
% 
% figure;
% imagesc(neiglinkmat);
% title(['Neighbor definition: ',neighbourdist ])
%  matcorr = zeros(83,83);
%                     matcorr(idhalf)=  neiglinkmat(:,329); %idhalf ;
%                  
%                     matcorr = matcorr +flipud(rot90(matcorr));
%                     meancorr = matcorr;

%%%% FIN to find neiglinkmat
%     neiglinkmat(21,:)=0
%     neiglinkmat(:,21)=0
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.connectivity     = neiglinkmat;
cfg.latency          = 'all';
cfg.frequency        = 'all';
cfg.numrandomization = 'all';
cfg.dim              = [size(statobs,1),1,1];
cfg.dimord           = 'chan_freq_time';
cfg.correctm         = 'cluster';
cfg.clusterstatistic = 'maxsum';
cfg.clusteralpha     = clustercritval;
cfg.clustercritval   = clustercritval;
cfg.minnbchan        = str2num(minnbchan);      %nb element significatif conjoints pour le cluster + haut car beaucoup de canaux
try
[stat, cfg] = clusterstat(cfg,statrand, statobs); 
catch 
    stat = []; 
    'no cluster'
end 
end 

function [neighbours] = compneighbstructfromgradelec(chanpos, label, neighbourdist)

nsensors = length(label);

% compute the distance between all sensors
dist = zeros(nsensors,nsensors);
for i=1:nsensors
  dist(i,:) = sqrt(sum((chanpos(1:nsensors,:) - repmat(chanpos(i,:), nsensors, 1)).^2,2))';
end

% find the neighbouring electrodes based on distance
% later we have to restrict the neighbouring electrodes to those actually selected in the dataset
channeighbstructmat = (dist<neighbourdist);

% electrode istelf is not a neighbour
channeighbstructmat = (channeighbstructmat .* ~eye(nsensors));

% construct a structured cell array with all neighbours
neighbours=struct;
for i=1:nsensors
  neighbours(i).label       = label{i};
  neighbours(i).neighblabel = label(find(channeighbstructmat(i,:)));
  neighbours(i).idneigh = find(channeighbstructmat(i,:));
end
end



function [sel1, sel2] = match_str(a, b, fullout)

% MATCH_STR looks for matching labels in two lists of strings
% and returns the indices into both the 1st and 2nd list of the matches.
% They will be ordered according to the first input argument.
%
% Use as
%   [sel1, sel2] = match_str(strlist1, strlist2)
%
% The strings can be stored as a char matrix or as an vertical array of
% cells, the matching is done for each row.
%
% When including a 1 as the third input argument, the output lists of
% indices will be expanded to the size of the largest input argument.
% Entries that occur only in one of the two inputs will correspond to a 0
% in the output, in this case. This can be convenient in rare cases if the
% size of the input lists is meaningful.

% Copyright (C) 2000-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% ensure that both are cell-arrays
if isempty(a)
  a = {};
elseif ~iscell(a)
  a = cellstr(a);
end
if isempty(b)
  b = {};
elseif ~iscell(b)
  b = cellstr(b);
end

% regardless of what optimizations are implemented, the code should remain
% functionally compatible to the original, which is
% sel1 = [];
% sel2 = [];
% for i=1:length(a)
%   for j=1:length(b)
%     if strcmp(a(i),b(j))
%        sel1 = [sel1; i];
%        sel2 = [sel2; j];
%      end
%    end
% end

% ensure that both are column vectors
a = a(:);
b = b(:);
Na = numel(a);
Nb = numel(b);

% According to the original implementation empty numeric elements are
% allowed, but are not returned as match. This is different to empty string
% elements, which are returned as match.
% See also http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1808
empty_a = cellfun(@isnumeric, a) & cellfun(@isempty, a);
empty_b = cellfun(@isnumeric, b) & cellfun(@isempty, b);
% the following allows the unique function to operate normally
a(empty_a) = {''};
b(empty_b) = {''};

% replace all unique strings by a unique number and use the fact that
% numeric comparisons are much faster than string comparisons
[dum1, dum2, c] = unique([a; b]);
a = c(1:Na);
b = c((Na+1):end);

% empty numeric elements should never be returned as a match
a(empty_a) = nan;
b(empty_b) = nan;

if nargin < 3 || ~fullout
  sel1 = [];
  sel2 = [];
  for i=1:length(a)
    % s = find(strcmp(a(i), b));  % for string comparison
    s = find(a(i)==b);            % for numeric comparison
    sel2 = [sel2; s];
    s(:) = i;
    sel1 = [sel1; s];
  end
else
  sel1 = zeros(max(Na,Nb),1);
  sel2 = zeros(max(Na,Nb),1);
  for i=1:length(a)
    s = find(a(i)==b);
    sel2(s) = s;
    sel1(s) = i;
  end
end
end

function [stat, cfg] = clusterstat(cfg, statrnd, statobs)

% CLUSTERSTAT computers cluster statistic for multidimensional channel-freq-time or
% volumetric source data
%
% See also TFCESTAT, FINDCLUSTER

% Copyright (C) 2005-2020, Robert Oostenveld
% Copyright (C) 2021, Robert Oostenveld and Jan-Mathijs Schoffelen
%  
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% set the defaults
cfg.feedback     = ft_getopt(cfg, 'feedback',     'text');
cfg.spmversion   = ft_getopt(cfg, 'spmversion',   'spm12');
cfg.dim          = ft_getopt(cfg, 'dim',          []);        % 1x3 vector, for volumetric source data
cfg.inside       = ft_getopt(cfg, 'inside',       []);
cfg.tail         = ft_getopt(cfg, 'tail',         0);         % -1, 0, 1

cfg.orderedstats = ft_getopt(cfg, 'orderedstats', 'no');      % no, yes
cfg.multivariate = ft_getopt(cfg, 'multivariate', 'no');      % no, yes
cfg.minnbchan    = ft_getopt(cfg, 'minnbchan',    0);
cfg.wcm_weight   = ft_getopt(cfg, 'wcm_weight',   1);

% these defaults are already set in the caller function, 
% but may be necessary if a user calls this function directly
cfg.clusterstatistic = ft_getopt(cfg, 'clusterstatistic', 'maxsum');      % max, maxsize, maxsum, wcm
cfg.clusterthreshold = ft_getopt(cfg, 'clusterthreshold', 'parametric');  % parametric, nonparametric_individual, nonparametric_common
cfg.clusteralpha     = ft_getopt(cfg, 'clusteralpha',     0.05);
cfg.clustercritval   = ft_getopt(cfg, 'clustercritval',   []);
cfg.clustertail      = ft_getopt(cfg, 'clustertail',      cfg.tail);
cfg.connectivity     = ft_getopt(cfg, 'connectivity',     false);

% ensure that the preferred SPM version is on the path
%ft_hastoolbox(cfg.spmversion, 1);

if isempty(cfg.dim)
  ft_error('cfg.dim should be defined and not empty');
end

if cfg.tail~=cfg.clustertail
  ft_error('cfg.tail and cfg.clustertail should be identical')
end

if isempty(cfg.inside)
  cfg.inside = true(cfg.dim);
end % cfg.inside is set in ft_sourcestatistics, but is also needed for timelock and freq

if isfield(cfg, 'origdim')
  cfg.dim = cfg.origdim;
end % this snippet is to support correct clustering of N-dimensional data, not fully tested yet

% get connectivity matrix for the spatially neighbouring elements
connmat = full(ft_getopt(cfg, 'connectivity', false));

needpos = cfg.tail==0 || cfg.tail== 1;
needneg = cfg.tail==0 || cfg.tail==-1;

Nsample    = size(statrnd,1);
Nrand      = size(statrnd,2);
prb_pos    = ones(Nsample,     1);
prb_neg    = ones(Nsample,     1);
postailrnd = false(Nsample,Nrand);  % this holds the thresholded values
negtailrnd = false(Nsample,Nrand);  % this holds the thresholded values
Nobspos    = 0;                     % number of positive clusters in observed data
Nobsneg    = 0;                     % number of negative clusters in observed data

switch cfg.clusterthreshold
  case 'parametric'
    if isempty(cfg.clustercritval)
      ft_error('with parametric cluster thresholding cfg.clustercritval needs to be defined');
    end
    
    % threshold based on the critical value from parametric distribution
    siz = size(cfg.clustercritval);
    if all(siz==1) && cfg.clustertail==0
      %  it only specifies one critical value, assume that the left and right tail are symmetric around zero
      negtailcritval = -cfg.clustercritval;
      postailcritval =  cfg.clustercritval;
    elseif all(siz==1) && cfg.clustertail==-1
      % it only specifies one critical value corresponding to the left tail
      negtailcritval = cfg.clustercritval;
      postailcritval = +inf * ones(size(negtailcritval));
    elseif all(siz==1) && cfg.clustertail==1
      % it only specifies one critical value corresponding to the right tail
      postailcritval =  cfg.clustercritval;
      negtailcritval = -inf * ones(size(postailcritval));
    elseif siz(1)==Nsample && siz(2)==1 && cfg.clustertail==0
      %  it specifies a single critical value for each sample, assume that the left and right tail are symmetric around zero
      negtailcritval = -cfg.clustercritval;
      postailcritval =  cfg.clustercritval;
    elseif siz(1)==Nsample && siz(2)==1 && cfg.clustertail==-1
      % it specifies a critical value for the left tail
      % which is different for each sample (samples have a different df)
      negtailcritval = cfg.clustercritval;
      postailcritval = +inf * ones(size(negtailcritval));
    elseif siz(1)==Nsample && siz(2)==1 && cfg.clustertail==1
      % it specifies a critical value for the right tail
      % which is different for each sample (samples have a different df)
      postailcritval = cfg.clustercritval;
      negtailcritval = +inf * ones(size(postailcritval));
    elseif siz(1)==Nsample && siz(2)==2 && cfg.clustertail==0
      % it specifies a critical value for the left and for the right tail of the distribution
      % which is different for each sample (samples have a different df)
      negtailcritval = cfg.clustercritval(:,1);
      postailcritval = cfg.clustercritval(:,2);
    elseif prod(siz)==2 && cfg.clustertail==0
      % it specifies a critical value for the left and for the right tail of the distribution
      % which is the same for each sample (samples have the same df)
      negtailcritval = cfg.clustercritval(1);
      postailcritval = cfg.clustercritval(2);
    else
      ft_error('cannot make sense out of the specified parametric critical values');
    end
    
  case 'nonparametric_individual'
    if isempty(cfg.clusteralpha)
      ft_error('with nonparametric_indivdual cluster thresholding cfg.clusteralpha needs to be defined');
    end
    
    % threshold based on bootstrap using all other randomizations
    % each voxel will get an individual threshold
    [srt, ind] = sort(statrnd,2);
    if cfg.clustertail==0
      % both tails are needed
      negtailcritval = srt(:,round((  cfg.clusteralpha/2)*size(statrnd,2)));
      postailcritval = srt(:,round((1-cfg.clusteralpha/2)*size(statrnd,2)));
    elseif cfg.clustertail==1
      % only positive tail is needed
      postailcritval = srt(:,round((1-cfg.clusteralpha)*size(statrnd,2)));
      negtailcritval = -inf * ones(size(postailcritval));
    elseif cfg.clustertail==-1
      % only negative tail is needed
      negtailcritval = srt(:,round((  cfg.clusteralpha)*size(statrnd,2)));
      postailcritval = +inf * ones(size(negtailcritval));
    end
    
  case 'nonparametric_common'
    if isempty(cfg.clusteralpha)
      ft_error('with nonparametric_common cluster thresholding cfg.clusteralpha needs to be defined');
    end
    
    % threshold based on bootstrap using all other randomizations
    % all voxels will get a common threshold
    srt = sort(statrnd(:));
    if cfg.clustertail==0
      % both tails are needed
      negtailcritval = srt(round((  cfg.clusteralpha/2)*numel(statrnd)));
      postailcritval = srt(round((1-cfg.clusteralpha/2)*numel(statrnd)));
    elseif cfg.clustertail==1
      % only positive tail is needed
      postailcritval = srt(round((1-cfg.clusteralpha)*numel(statrnd)));
      negtailcritval = -inf * ones(size(postailcritval));
    elseif cfg.clustertail==-1
      % only negative tail is needed
      negtailcritval = srt(round((  cfg.clusteralpha)*numel(statrnd)));
      postailcritval = +inf * ones(size(negtailcritval));
    end
  case 'none'
    % tfce
    negtailcritval = [];
    postailcritval = [];
  otherwise
    ft_error('no valid threshold for clustering was given')
end % determine clusterthreshold

% these should be scalars or column vectors
negtailcritval = negtailcritval(:);
postailcritval = postailcritval(:);

% remember the critical values
cfg.clustercritval = [negtailcritval postailcritval];

% test whether the observed and the random statistics exceed the threshold
postailobs = (statobs >= postailcritval);
negtailobs = (statobs <= negtailcritval);
for i=1:Nrand
  postailrnd(:,i) = (statrnd(:,i) >= postailcritval);
  negtailrnd(:,i) = (statrnd(:,i) <= negtailcritval);
end

% first do the clustering on the observed data
spacereshapeable = (numel(connmat)==1 && ~isfinite(connmat));

if needpos
  if spacereshapeable
    % this pertains to data for which the spatial dimension can be reshaped
    % into 3D, i.e. when it is described on an ordered set of positions on
    % a 3D-grid. It deals with the inside dipole positions, and creates a
    % fake extra spatial dimension, so that findcluster can deal with it
    tmp = zeros([1 cfg.dim]); 
    tmp(cfg.inside) = postailobs; 
  else
    tmp = reshape(postailobs, [cfg.dim 1]);
  end
  
  % identify positive clusters in the observed data
  posclusobs = findcluster(tmp, connmat, cfg.minnbchan);
  
  if spacereshapeable
    posclusobs = posclusobs(cfg.inside);
  else
    posclusobs = posclusobs(:);
  end
  Nobspos    = max(posclusobs); % number of clusters exceeding the threshold
  fprintf('found %d positive clusters in observed data\n', Nobspos);
  
end % if needpos

if needneg
  if spacereshapeable
    % this pertains to data for which the spatial dimension can be reshaped
    % into 3D, i.e. when it is described on an ordered set of positions on
    % a 3D-grid. It deals with the inside dipole positions, and creates a
    % fake extra spatial dimension, so that findcluster can deal with it
    tmp = zeros([1 cfg.dim]); 
    tmp(cfg.inside) = negtailobs; 
  else
    tmp = reshape(negtailobs, [cfg.dim 1]);
  end
  
  % identify negative clusters in the observed data
  negclusobs = findcluster(tmp, connmat, cfg.minnbchan);
  
  if spacereshapeable
    negclusobs = negclusobs(cfg.inside);
  else
    negclusobs = negclusobs(:);
  end
  Nobsneg    = max(negclusobs); % number of clusters exceeding the threshold
  fprintf('found %d negative clusters in observed data\n', Nobsneg);
  
end % if needneg

% catch situation where no clustering of the random data is needed
if (Nobspos+Nobsneg)==0
  ft_warning('no clusters were found in the observed data');
  stat = struct(); % see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2972
  stat.stat = statobs;
  stat.prob = ones(Nsample, 1);
  return
end

% allocate space to hold the randomization distributions of the cluster statistic
if strcmp(cfg.multivariate, 'yes') || strcmp(cfg.orderedstats, 'yes')
  fprintf('allocating space for a %d-multivariate distribution of the positive clusters\n', Nobspos);
  fprintf('allocating space for a %d-multivariate distribution of the negative clusters\n', Nobsneg);
  posdistribution = zeros(Nobspos,Nrand); % this holds the multivariate randomization distribution of the positive cluster statistics
  negdistribution = zeros(Nobsneg,Nrand); % this holds the multivariate randomization distribution of the negative cluster statistics
else
  posdistribution = zeros(1,Nrand);       % this holds the statistic of the largest positive cluster in each randomization
  negdistribution = zeros(1,Nrand);       % this holds the statistic of the largest negative cluster in each randomization
end

% do the clustering on the randomized data
ft_progress('init', cfg.feedback, 'computing clusters for the thresholded test statistic computed from the randomized design');
for i = 1:Nrand
  ft_progress(i/Nrand, 'computing clusters in randomization %d from %d\n', i, Nrand);
  if needpos
    if spacereshapeable
      tmp = zeros([1 cfg.dim]);
      tmp(cfg.inside) = postailrnd(:,i);
    else
      tmp = reshape(postailrnd(:,i), [cfg.dim 1]);
    end
    posclusrnd = findcluster(tmp, connmat, cfg.minnbchan);
    if spacereshapeable
      posclusrnd = posclusrnd(cfg.inside);
    else
      posclusrnd = posclusrnd(:);
    end
    Nrndpos = max(posclusrnd(:));  % number of clusters exceeding the threshold
    stat    = zeros(1,Nrndpos); % this will hold the statistic for each cluster
    % fprintf('found %d positive clusters in this randomization\n', Nrndpos);
    for j = 1:Nrndpos
      switch cfg.clusterstatistic
        case 'max'
          stat(j) = max(statrnd(posclusrnd==j,i));
        case 'maxsize'
          stat(j) = length(find(posclusrnd==j));
        case 'maxsum'
          stat(j) = sum(statrnd(posclusrnd==j,i));
        case 'wcm'
          if numel(postailcritval)==1
            posthr = postailcritval;
          elseif numel(postailcritval)==numel(posclusrnd)
            posthr = postailcritval(posclusrnd==j);
          end
          stat(j) = sum((statrnd(posclusrnd==j,i)-posthr).^cfg.wcm_weight);
        otherwise
          ft_error('unknown clusterstatistic');
      end
    end % for 1:Nrdnpos
    if strcmp(cfg.multivariate, 'yes') || strcmp(cfg.orderedstats, 'yes')
      stat = sort(stat, 'descend');             % sort them from most positive to most negative
      if Nrndpos>Nobspos
        posdistribution(:,i) = stat(1:Nobspos); % remember the largest N clusters
      else
        posdistribution(1:Nrndpos,i) = stat;    % remember the largest N clusters
      end
    else
      % univariate -> remember the most extreme cluster
      if ~isempty(stat), posdistribution(i) = max(stat); end
    end
  end % needpos
  if needneg
    if spacereshapeable
      tmp = zeros([1 cfg.dim]);
      tmp(cfg.inside) = negtailrnd(:,i);
    else
      tmp = reshape(negtailrnd(:,i), [cfg.dim 1]);
    end
    negclusrnd = findcluster(tmp, connmat, cfg.minnbchan);
    if spacereshapeable
      negclusrnd = negclusrnd(cfg.inside);
    else  
      negclusrnd = negclusrnd(:);
    end
    Nrndneg = max(negclusrnd(:));  % number of clusters exceeding the threshold
    stat    = zeros(1,Nrndneg); % this will hold the statistic for each cluster
    % fprintf('found %d negative clusters in this randomization\n', Nrndneg);
    for j = 1:Nrndneg
      switch cfg.clusterstatistic
        case 'max'
          stat(j) = min(statrnd(negclusrnd==j,i));
        case 'maxsize'
          stat(j) = -length(find(negclusrnd==j)); % encode the size of a negative cluster as a negative value
        case 'maxsum'
          stat(j) = sum(statrnd(negclusrnd==j,i));
        case 'wcm'
          if numel(negtailcritval)==1
            negthr = negtailcritval;
          elseif numel(negtailcritval)==numel(negclusrnd)
            negthr = negtailcritval(negclusrnd==j);
          end
          stat(j) = -sum((abs(statrnd(negclusrnd==j,i)-negthr)).^cfg.wcm_weight); % encoded as a negative value
        otherwise
          ft_error('unknown clusterstatistic');
      end
    end % for 1:Nrndneg
    if strcmp(cfg.multivariate, 'yes') || strcmp(cfg.orderedstats, 'yes')
      stat = sort(stat, 'ascend');              % sort them from most negative to most positive
      if Nrndneg>Nobsneg
        negdistribution(:,i) = stat(1:Nobsneg); % remember the most extreme clusters, i.e. the most negative
      else
        negdistribution(1:Nrndneg,i) = stat;    % remember the most extreme clusters, i.e. the most negative
      end
    else
      % univariate -> remember the most extreme cluster, which is the most negative
      if ~isempty(stat), negdistribution(i) = min(stat); end
    end
  end % needneg
end % for 1:Nrand
ft_progress('close');

% compare the values for the observed clusters with the randomization distribution
if needpos
  posclusters = [];
  stat = zeros(1,Nobspos);
  for j = 1:Nobspos
    switch cfg.clusterstatistic
      case 'max'
        stat(j) = max(statobs(posclusobs==j));
      case 'maxsize'
        stat(j) = length(find(posclusobs==j));
      case 'maxsum'
        stat(j) = sum(statobs(posclusobs==j));
      case 'wcm'
        if numel(postailcritval)==1
          posthr = postailcritval;
        elseif numel(postailcritval)==numel(posclusrnd)
          posthr = postailcritval(posclusobs==j);
        end
        stat(j) = sum((statobs(posclusobs==j)-posthr).^cfg.wcm_weight);
      otherwise
        ft_error('unknown clusterstatistic');
    end
  end
  % sort the clusters based on their statistical value
  [stat, indx] = sort(stat, 'descend');
  % reorder the cluster indices in the data
  tmp = zeros(size(posclusobs));
  for j=1:Nobspos
    tmp(posclusobs==indx(j)) = j;
  end
  posclusobs = tmp;
  if strcmp(cfg.multivariate, 'yes')
    % estimate the probability of the mutivariate tail, i.e. one p-value for all clusters
    prob = 0;
    for i=1:Nrand
      % compare all clusters simultaneosuly
      prob = prob + any(posdistribution(:,i)>stat(:));
    end
    if isequal(cfg.numrandomization, 'all')
      prob = prob/Nrand;
    else % the minimum possible p-value should not be 0, but 1/N
      prob = (prob + 1)/(Nrand + 1);
    end
    for j = 1:Nobspos
      % collect a summary of the cluster properties
      posclusters(j).prob = prob;
      posclusters(j).clusterstat = stat(j);
    end
    % collect the probabilities in one large array
    prb_pos(posclusobs~=0) = prob;
  elseif strcmp(cfg.orderedstats, 'yes')
    % compare the Nth ovbserved cluster against the randomization distribution of the Nth cluster
    prob = zeros(1,Nobspos);
    for j = 1:Nobspos
      if isequal(cfg.numrandomization, 'all')
        prob(j) = sum(posdistribution(j,:)>stat(j))/Nrand;
      else % the minimum possible p-value should not be 0, but 1/N
        prob(j) = (sum(posdistribution(j,:)>stat(j)) + 1)/(Nrand + 1);
      end
      % collect the probabilities in one large array
      prb_pos(posclusobs==j) = prob(j);
    end
    % collect a summary of the cluster properties
    posclusters = struct('prob', num2cell(prob), 'clusterstat', num2cell(stat));    
  else
    % univariate -> each cluster has it's own probability
    prob = zeros(1,Nobspos);
    for j = 1:Nobspos
      if isequal(cfg.numrandomization, 'all')
        prob(j) = sum(posdistribution>stat(j))/Nrand;
      else % the minimum possible p-value should not be 0, but 1/N
        prob(j) = (sum(posdistribution>stat(j)) + 1)/(Nrand + 1);
      end
      % collect the probabilities in one large array
      prb_pos(posclusobs==j) = prob(j);
    end
    % collect a summary of the cluster properties
    posclusters = struct('prob', num2cell(prob), 'clusterstat', num2cell(stat));    
  end
end

if needneg
  negclusters = [];
  stat = zeros(1,Nobsneg);
  for j = 1:Nobsneg
    switch cfg.clusterstatistic
      case 'max'
        stat(j) = min(statobs(negclusobs==j));
      case 'maxsize'
        stat(j) = -length(find(negclusobs==j)); % encode the size of a negative cluster as a negative value
      case 'maxsum'
        stat(j) = sum(statobs(negclusobs==j));
      case 'wcm'
        if numel(negtailcritval)==1
          negthr = negtailcritval;
        elseif numel(negtailcritval)==numel(negclusrnd)
          negthr = negtailcritval(negclusobs==j);
        end
        stat(j) = -sum((abs(statobs(negclusobs==j)-negthr)).^cfg.wcm_weight); % encoded as a negative value
      otherwise
        ft_error('unknown clusterstatistic');
    end
  end
  % sort the clusters based on their statistical value
  [stat, indx] = sort(stat,'ascend');
  % reorder the cluster indices in the observed data
  tmp = zeros(size(negclusobs));
  for j=1:Nobsneg
    tmp(negclusobs==indx(j)) = j;
  end
  negclusobs = tmp;
  if strcmp(cfg.multivariate, 'yes')
    % estimate the probability of the mutivariate tail, i.e. one p-value for all clusters
    prob = 0;
    for i=1:Nrand
      % compare all clusters simultaneosuly
      prob = prob + any(negdistribution(:,i)<stat(:));
    end
    if isequal(cfg.numrandomization, 'all')
      prob = prob/Nrand;
    else % the minimum possible p-value should not be 0, but 1/N
      prob = (prob + 1)/(Nrand + 1);
    end
    for j = 1:Nobsneg
      % collect a summary of the cluster properties
      negclusters(j).prob = prob;
      negclusters(j).clusterstat = stat(j);
    end
    % collect the probabilities in one large array
    prb_neg(negclusobs~=0) = prob;
  elseif strcmp(cfg.orderedstats, 'yes')
    % compare the Nth ovbserved cluster against the randomization distribution of the Nth cluster
    prob = zeros(1,Nobsneg);
    for j = 1:Nobsneg
      if isequal(cfg.numrandomization, 'all')
        prob(j) = sum(negdistribution(j,:)<stat(j))/Nrand;
      else % the minimum possible p-value should not be 0, but 1/N
        prob(j) = (sum(negdistribution(j,:)<stat(j)) + 1)/(Nrand + 1);
      end
      % collect the probabilities in one large array
      prb_neg(negclusobs==j) = prob(j);
    end
    % collect a summary of the cluster properties
    negclusters = struct('prob', num2cell(prob), 'clusterstat', num2cell(stat));
  else
    % univariate -> each cluster has its own probability
    prob = zeros(1,Nobsneg);
    for j = 1:Nobsneg
      if isequal(cfg.numrandomization, 'all')
        prob(j) = sum(negdistribution<stat(j))/Nrand;
      else % the minimum possible p-value should not be 0, but 1/N
        prob(j) = (sum(negdistribution<stat(j)) + 1)/(Nrand + 1);
      end
      % collect the probabilities in one large array
      prb_neg(negclusobs==j) = prob(j);
    end
    % collect a summary of the cluster properties
    negclusters = struct('prob', num2cell(prob), 'clusterstat', num2cell(stat));
  end
end

if cfg.tail==0
  % consider both tails
  prob = min(prb_neg, prb_pos); % this is the probability for the most unlikely tail
elseif cfg.tail==1
  % only consider the positive tail
  prob = prb_pos;
elseif cfg.tail==-1
  % only consider the negative tail
  prob = prb_neg;
end

% collect the remaining details in the output structure
stat = struct(); % see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2972
stat.prob = prob;
if needpos
  stat.posclusters         = posclusters;
  stat.posclusterslabelmat = posclusobs;
  stat.posdistribution     = posdistribution;
end
if needneg
  stat.negclusters         = negclusters;
  stat.negclusterslabelmat = negclusobs;
  stat.negdistribution     = negdistribution;
end

end
function [cluster, numcluster] = findcluster(onoff, spatdimneighbstructmat, varargin)

% FINDCLUSTER returns all connected clusters for a three-dimensional six-connected
% neighborhood
%
% Use as
%   [cluster, num] = findcluster(onoff, spatdimneighbstructmat, minnbchan)
% or ar
%   [cluster, num] = findcluster(onoff, spatdimneighbstructmat, spatdimneighbselmat, minnbchan)
% where
%   onoff                  =  is a 3D boolean matrix with size N1xN2xN3
%   spatdimneighbstructmat =  defines the neighbouring channels/combinations, see below
%   minnbchan              =  the minimum number of neighbouring channels/combinations
%   spatdimneighbselmat    =  is a special neighbourhood matrix that is used for selecting
%                             channels/combinations on the basis of the minnbchan criterium
%
% The neighbourhood structure for the first dimension is specified using
% spatdimneighbstructmat, which is a 2D (N1xN1) matrix. Each row and each column
% corresponds to a channel (combination) along the first dimension and along that
% row/column, elements with "1" define the neighbouring channel(s) (combinations).
% The first dimension of onoff should correspond to the channel(s) (combinations).
%
% See also SPM_BWLABEL, BWLABEL, BWLABELN

% Copyright (C) 2004-2020, Robert Oostenveld
% Copyright (C) 2021, Robert Oostenveld and Jan-Mathijs Schoffelen
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% the calling code should ensure that SPM is on the path, preferably the latest version
%ft_hastoolbox('spm', -1);

siz           = size(onoff);
spatdimlength = siz(1);
siz           = siz(2:end);
dims          = ndims(onoff);

% this is for reshaping back
celldims      = num2cell([siz 1]);

% this is for efficient matrix indexing
seldims = repmat({':'},[1 dims-1]);

% put the spatial dimension last: this substantially speeds up allocation
% and selection of data from the matrix, which outweighs the time spent to
% do the matrix permutation;
onoff = permute(onoff, [2:dims 1]);

if ~ismatrix(spatdimneighbstructmat) || ~all(size(spatdimneighbstructmat)==spatdimlength)
  ft_error('invalid dimension of spatdimneighbstructmat');
end

% this input argument handling is somewhat clunky, but it does the trick
if length(varargin)==1
  minnbchan = varargin{1};
else
  minnbchan = 0;
end
if length(varargin)==2
  spatdimneighbselmat = varargin{1};
  minnbchan           = varargin{2};
end

if minnbchan>0
  % For every (time,frequency)-element, it is calculated how many supra
  % threshold neighbours this spatial element has. If a supra threshold
  % spatial elementhas fewer than minnbchan supra threshold neighbours,
  % it is removed from onoff.
  
  if length(varargin)==1
    selectmat = single(spatdimneighbstructmat | spatdimneighbstructmat');
  end
  if length(varargin)==2
    selectmat = single(spatdimneighbselmat | spatdimneighbselmat');
  end
  
  nremoved = 1;
  while nremoved>0
    nsigneighb = reshape(reshape(single(onoff),[prod(siz) spatdimlength])*selectmat,[siz spatdimlength]);
    remove     = (onoff.*nsigneighb) < minnbchan;
    nremoved   = length(find(remove.*onoff));
    onoff(remove) = 0;
  end
end

% for each channel or channel-combination, find the connected time, frequency, or time-frequency clusters
labelmat = zeros(size(onoff));
numcluster = 0;

if ~(numel(siz)==1 && all(siz==1) && islogical(onoff))
  for j = 1:spatdimlength
    if numel(siz) <= 3 % if 2D or 3D data (without spatial dimension)
      % use SPM for 2D/3D data instead of the MATLAB image processing toolbox
      [clus, num] = spm_bwlabel(double(onoff(seldims{:},j)), 6);
    else
      [clus, num] = bwlabeln(double(onoff(seldims{:},j)), conndef(dims-1, 'min'));
    end
    clus(clus~=0) = clus(clus~=0) + numcluster;
    labelmat(seldims{:},j) = clus;
    numcluster = numcluster + num;
  end
else
  labelmat(onoff>0) = 1:sum(onoff(:));
  numcluster = sum(onoff(:));
end

% put the spatial dimension back upfront
labelmat = permute(labelmat, [dims 1:(dims-1)]);

% combine the non-spatial dimensions for simplicity
labelmat = reshape(labelmat, spatdimlength, []);

% combine clusters that are connected in neighbouring channels or channel
% combinations. Here we convert the input to uint32 as that is required by the mex
% file, and the values will be positive integers anyway.
if spatdimlength>1
  cluster = combineClusters(uint32(labelmat), logical(spatdimneighbstructmat), uint32(numcluster));
else
  cluster = labelmat;
end

% reshape the output to the original format of the data
cluster = reshape(cluster, spatdimlength, celldims{:});

% update the total number of clusters
numcluster = max(cluster(:));
end
function cluster = combineClusters(labelmat, spatdimneighbstructmat, total)

% COMBINECLUSTERS is a helper function for FINDCLUSTER. It searches for
% adjacent clusters in neighbouring channels and combines them.

% A mex-file is available for this function, as it can take quite long.

replaceby=1:total;
spatdimlength = size(labelmat, 1);
for spatdimlev=1:spatdimlength
  neighbours=find(spatdimneighbstructmat(spatdimlev,:));
  for nbindx=neighbours
    indx = find((labelmat(spatdimlev,:)~=0) & (labelmat(nbindx,:)~=0));
    for i=1:length(indx)
      a = labelmat(spatdimlev, indx(i));
      b = labelmat(nbindx, indx(i));
      if replaceby(a)==replaceby(b)
        % do nothing
        continue;
      elseif replaceby(a)<replaceby(b)
        % replace all entries with content replaceby(b) by replaceby(a).
        replaceby(find(replaceby==replaceby(b))) = replaceby(a); 
      elseif replaceby(b)<replaceby(a)
        % replace all entries with content replaceby(a) by replaceby(b).
        replaceby(find(replaceby==replaceby(a))) = replaceby(b); 
      end
    end
  end
end

% renumber all the clusters
num = 0;
cluster = zeros(size(labelmat));
for uniquelabel=unique(replaceby(:))'
  num = num+1;
  cluster(ismember(labelmat(:),find(replaceby==uniquelabel))) = num;
end

end

function val = ft_getopt(opt, key, default, emptymeaningful)

% FT_GETOPT gets the value of a specified option from a configuration structure
% or from a cell-array with key-value pairs.
%
% Use as
%   val = ft_getopt(s, key, default, emptymeaningful)
% where the input values are
%   s               = structure or cell-array
%   key             = string
%   default         = any valid MATLAB data type (optional, default = [])
%   emptymeaningful = boolean value (optional, default = 0)
%
% If the key is present as field in the structure, or as key-value pair in the
% cell-array, the corresponding value will be returned.
%
% If the key is not present, ft_getopt will return the default, or an empty array
% when no default was specified.
%
% If the key is present but has an empty value, then the emptymeaningful flag
% specifies whether the empty value or the default value should be returned.
% If emptymeaningful==true, then the empty array will be returned.
% If emptymeaningful==false, then the specified default will be returned.
%
% See also FT_SETOPT, FT_CHECKOPT

% Copyright (C) 2011-2012, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

if nargin<3
  default = [];
end

if nargin < 4
  emptymeaningful = 0;
end

if isa(opt, 'struct') || isa(opt, 'config')
  % get the key-value from the structure
  fn = fieldnames(opt);
  if ~any(strcmp(key, fn))
    val = default;
  else
    val = opt.(key);
  end
  
elseif isa(opt, 'cell')
  % get the key-value from the cell-array
  if mod(length(opt),2)
    error('optional input arguments should come in key-value pairs, i.e. there should be an even number');
  end
  
  % the 1st, 3rd, etc. contain the keys, the 2nd, 4th, etc. contain the values
  keys = opt(1:2:end);
  vals = opt(2:2:end);
  
  % the following may be faster than cellfun(@ischar, keys)
  valid = false(size(keys));
  for i=1:numel(keys)
    valid(i) = ischar(keys{i});
  end
  
  if ~all(valid)
    error('optional input arguments should come in key-value pairs, the optional input argument %d is invalid (should be a string)', i);
  end
  
  hit = find(strcmpi(key, keys));
  if isempty(hit)
    % the requested key was not found
    val = default;
  elseif length(hit)==1
    % the requested key was found
    val = vals{hit};
  else
    error('multiple input arguments with the same name');
  end
  
elseif isempty(opt)
  % no options are specified, return default
  val = default;
end % isstruct or iscell or isempty

if isempty(val) && ~isempty(default) && ~emptymeaningful
  % use the default value instead of the empty input that was specified:
  % this applies for example if you do functionname('key', []), where
  % the empty is meant to indicate that the user does not know or care
  % what the value is
  val = default;
end
end
function ft_progress(varargin)

% FT_PROGRESS shows a graphical or non-graphical progress indication similar to the
% standard WAITBAR function, but with the extra option of printing it in the command
% window as a plain text string or as a rotating dial. Alternatively, you can also
% specify it not to give feedback on the progress.
%
% Prior to the for-loop, you should call either
%   ft_progress('init', 'none',    'Please wait...')
%   ft_progress('init', 'gui',     'Please wait...')
%   ft_progress('init', 'etf',     'Please wait...')      % estimated time to finish
%   ft_progress('init', 'dial',    'Please wait...')      % rotating dial
%   ft_progress('init', 'textbar', 'Please wait...')      % ascii progress bar
%   ft_progress('init', 'text',    'Please wait...')
%
% In each iteration of the for-loop, you should call either
%   ft_progress(x)                                       % only show percentage
%   ft_progress(x, 'Processing event %d from %d', i, N)  % show string, x=i/N
%
% After finishing the for-loop, you should call
%   ft_progress('close')
%
% Here is an example for the use of a progress indicator
%    ft_progress('init', 'etf',     'Please wait...');
%    for i=1:42
%      ft_progress(i/42, 'Processing event %d from %d', i, 42);
%      pause(0.1);
%    end
%    ft_progress('close')
%
% See also WAITBAR

% Copyright (C) 2004-2008, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

global ft_default

persistent p        % the previous value of the progress
persistent c        % counter for the number of updates that is done
persistent t0       % initial time, required for ETF
persistent p0       % initial percentage, required for ETF
persistent t        % type of feedback, string with none, gui, text, textcr, textnl
persistent h        % the handle of the dialog (in case of type=gui)
persistent a        % the angle in degrees, for dial or textbar
persistent s        % the string containing the title
persistent strlen   % the length of the previously printed string, used to remove it by \b
persistent tprev    % the time of previous invocation, used to restrict number of updates

persistent lastArgin % the last varargin, this is used when ft_progress('close') is called
                     % but the previous invocation was not processed (due
                     % to the restriction in the number of updates to once
                     % every 100ms)
persistent closing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin>1 && ischar(varargin{1}) && strcmp(varargin{1}, 'init')
  % reset these to the defaults
  a  = 0;
  h  = 0;
  p  = 0;
  t  = 'none';
  s  = '';
  t0 = [];
  p0 = [];
  strlen = 0;
  tprev = tic();
  lastArgin = [];
  closing = 0;
  ft_default.progress.noerase = 0;
  
  % determine the type of feedback
  t = varargin{2};
  if strcmp(t,'textcr') || strcmp(t,'textnl')
    t = 'text';
  end
  % determine the title of the dialog
  if nargin>2
    s = varargin{3};
  else
    s = '';
  end
  switch t
  case 'gui'
    % initialise the waitbar dialog
    if ~isempty(s)
      h = waitbar(0, s);
    else
      h = waitbar(0, 'Please wait');
    end
  case {'text', 'textnl', 'textcr', 'dial'}
    if ~isempty(s)
      % print the title to the screen and go to the next line
      fprintf('%s\n', s)
    end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif nargin==1 && ischar(varargin{1}) && strcmp(varargin{1}, 'close')
  
  if ~isempty(lastArgin)
    closing = 1;
    ft_progress(lastArgin{:});
  end
  
  switch t
  case 'gui'
    % close the waitbar dialog
    close(h);
  case {'text', 'etf', 'dial', 'textbar'}
    % finish by going to the next line
    fprintf('\n');
  end
  % reset these to the defaults
  a  = 0;
  h  = 0;
  p  = 0;
  t  = 'none';
  s  = '';
  t0 = [];
  p0 = [];
  strlen = 0;
  tprev = [];
  lastArgin = [];
  closing = 0;
  ft_default.progress.noerase = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  
  if isfield(ft_default, 'progress') &&...
      isfield(ft_default.progress, 'noerase') &&...
      ft_default.progress.noerase
    strlen = 0;
    ft_default.progress.noerase = 0;
  end
  
  % make sure we don't update more than once every 100ms, significant
  % performance hit otherwise in certain conditions
  if ~isempty(tprev) && toc(tprev) < 0.1 && ~closing
    lastArgin = varargin;
    return;
  end
  tprev = tic();
  lastArgin = [];
  
  if strcmp(t, 'dial')
    % display should always be updated for the dial
    % continue;
  elseif (varargin{1}-p)<0.01 && strcmp(t, 'gui')
    % display should not be updated it the difference is less than one percent
    return;
  elseif (varargin{1}-p)<0.01 && strcmp(t, 'textbar')
    % display should not be updated it the difference is less than one percent
    return;
  elseif (varargin{1}-p)<0.01 && strcmp(t, 'etf')
    % display should not be updated it the difference is less than one percent
    % return;
  end

  % count the number of updates, for debugging
  c = c+1;

  % remember the current value for the next function call
  p = varargin{1};

  switch t
  case 'gui'
    % update the the length of the bar in the waitbar dialog
    waitbar(varargin{1}, h);

  case 'etf'
    % compute the estimated time that the computation still needs to finish
    if isempty(t0) || isempty(p0)
      t0 = clock;
      p0 = p;
    end
    elapsed = etime(clock, t0);
    if nargin>1 && ~isempty(varargin{2})
      
      varargin{2} = [repmat(sprintf('\b'),[1 strlen]) varargin{2}];
      
      % include the specified string
      strlentmp = fprintf(varargin{2:end});
      strlentmp = strlentmp + fprintf(' - estimated time to finish is %d seconds', round(elapsed*(1-p)/(p-p0)));
      
      % record actual string length that was printed (subtracting all the
      % \b's)
      strlen = strlentmp - strlen;
    else
      % only print the estimated time to finish
      strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) ' - estimated time to finish is %d seconds'], round(elapsed*(1-p)/(p-p0)));
      strlen = strlentmp - strlen;
    end

  case 'dial'
    dial = '|/-\|/-\';
    if ~isempty(s)
      % print the title and draw a new hand of the rotating dial
      strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s %s'], s, dial(1+a/45));
      strlen = strlentmp - strlen;
    else
      % draw a new hand of the rotating dial
      strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%s'], dial(1+a/45));
      srtlen = strlentmp - strlen;
    end
    % increment the angle with 45 degrees
    a = a + 45;
    if a==360
      % reset the angle to 0 degrees
      a = 0;
    end

  case 'textbar'
    dial = '|/-\|/-\';
    % construct the line looking like [------/          ]
    len  = 75 - length(s) - 3;
    len1 = round(p*len);            % number of '-' characters before the dial
    len2 = len - len1;              % number of ' ' characters after the dial
    line = [s, ' [' repmat('-',1,len1), dial(1+a/45), repmat(' ',1,len2) ,']'];
    if c~=1
      backline = repmat('\b', [1 length(line)]);
    else
      backline = '';
    end
    fprintf([backline,'%s'], line);
    %fprintf('\r%s', line); %carriage return sometimes leads to a new line
    % increment the angle with 45 degrees
    a = a + 45;
    if a==360
      % reset the angle to 0 degrees
      a = 0;
    end

  case 'text'
    if nargin>1
      
      % ensure the string does not end with a newline or carriage return
      % either would break compatibility with a -nodesktop matlab
      % environment
      if length(varargin{2})>1 && (all(varargin{2}((end-1):end) == '\r')...
          || all(varargin{2}((end-1):end) == '\n'))
        varargin{2} = varargin{2}(1:end-2);
      end
      
      varargin{2} = [repmat(sprintf('\b'),[1 strlen]) varargin{2}];
      if usejava('desktop')
        % a newline is appropriate when using the desktop environment
        varargin{2} = [varargin{2} '\n'];
      end
      
      strlentmp = fprintf(varargin{2:end});
      strlen = strlentmp - strlen;
    else
      strlentmp = fprintf([repmat(sprintf('\b'),[1 strlen]) '%6.2f %%'], 100*varargin{1});
      strlen = strlentmp - strlen;
    end

% the following options are unused in FieldTrip (as of April 17 2012), and seem
% semantically incompatible with the implementation of the \b-ing, so I
% think removal is appropriate.
%
%   case 'textnl'
%     if nargin>1
%       % ensure that the string ends with a newline
%       if length(varargin{2})>1 && all(varargin{2}((end-1):end) == '\r')
%         varargin{2}((end-1):end) = '\n';
%       elseif length(varargin{2})>1 && ~all(varargin{2}((end-1):end) == '\n')
%         varargin{2}((end+1):(end+2)) = '\n';
%       elseif length(varargin{2})<2
%         varargin{2}((end+1):(end+2)) = '\n';
%       end
%       fprintf(varargin{2:end});
%     else
%       fprintf('%6.2f %%\n', 100*varargin{1});
%     end
% 
%   case 'textcr'
%     if nargin>1
%       % ensure that the string ends with a cariage return
%       if length(varargin{2})>1 && all(varargin{2}((end-1):end) == '\n')
%         varargin{2}((end-1):end) = '\r';
%       elseif length(varargin{2})>1 && ~all(varargin{2}((end-1):end) == '\r')
%         varargin{2}((end+1):(end+2)) = '\r';
%       elseif length(varargin{2})<2
%         varargin{2}((end+1):(end+2)) = '\r';
%       end
%       fprintf(varargin{2:end});
%     else
%       fprintf('%6.2f %%\r', 100*varargin{1});
%     end

  end % case gui, dial, text
end % updating the displayed value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some test code follows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% progress('init', 'gui');  for i=1:100, progress(i/100, '%d/%d\r', i, 100); end; progress('close')
% progress('init', 'dial'); for i=1:100, progress(i/100, '%d/%d\r', i, 100); end; progress('close')
% progress('init', 'none'); for i=1:100, progress(i/100, '%d/%d\r', i, 100); end; progress('close')

% progress('init', 'text');   for i=1:100, progress(i/100, '%d/%d'  , i, 100); end; progress('close')
% progress('init', 'text');   for i=1:100, progress(i/100, '%d/%d\n', i, 100); end; progress('close')
% progress('init', 'text');   for i=1:100, progress(i/100, '%d/%d\r', i, 100); end; progress('close')

% progress('init', 'textnl'); for i=1:100, progress(i/100, '%d/%d'  , i, 100); end; progress('close')
% progress('init', 'textnl'); for i=1:100, progress(i/100, '%d/%d\n', i, 100); end; progress('close')
% progress('init', 'textnl'); for i=1:100, progress(i/100, '%d/%d\r', i, 100); end; progress('close')

% progress('init', 'textcr'); for i=1:100, progress(i/100, '%d/%d'  , i, 100); end; progress('close')
% progress('init', 'textcr'); for i=1:100, progress(i/100, '%d/%d\n', i, 100); end; progress('close')
% progress('init', 'textcr'); for i=1:100, progress(i/100, '%d/%d\r', i, 100); end; progress('close')


end