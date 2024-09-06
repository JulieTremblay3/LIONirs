function srx2seg(file,nb_x,nb_y,nb_z)
% Use a brainsuite surface file to create an imagic segmentation file 
% fileSRX
% New file create file fileSRX.seg, fileSRX.hdr
SRXTMP = readdfs(file)
SRX.vertices(:,1) = nb_x-SRXTMP.vertices(:,1);
SRX.vertices(:,2) = nb_y-SRXTMP.vertices(:,2);
SRX.vertices(:,3) = SRXTMP.vertices(:,3);
volumeIRM = zeros(nb_x,nb_y,nb_z);
for j = round(max(SRX.vertices(:,2))):-1:(round(min(SRX.vertices(:,2)))+1) %Check de haut en bas
    list = find(SRX.vertices(:,2)>(j-1) & SRX.vertices(:,2)<j );
      %Contours     
     x = SRX.vertices(list,1);
     idx = round(min(x)):1:round(max(x));
     for i=round(min(x)):1:round(max(x));
        listk = find(SRX.vertices(list,1)>(i-1) & SRX.vertices(list,1)<i);        
        minz = round(min(SRX.vertices(list(listk),3)));
        if minz>nb_z/2
            minz=1;
        end            
        maxz= round(max(SRX.vertices(list(listk),3)));
        volumeIRM(i,j,minz:maxz)=255; 
     end
  end 

 [pathstr, name, ext]=fileparts(file);
 %SAVE HDR
 filehdr = fullfile(pathstr,[name,'.hdr']);
 fid=fopen(filehdr, 'w');
 fprintf(fid,'%d\r\n',nb_x);
 fprintf(fid,'%d\r\n',nb_z);
 fprintf(fid,'%d\r\n',nb_y);
 for i=0:nb_z-1
      fprintf(fid,'%d\r\n',i);
 end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%  
fileseg = fullfile(pathstr,[name,'.seg']);

fid = fopen(fileseg,'wb');
fwrite(fid,volumeIRM,'uint8');
fclose(fid);

