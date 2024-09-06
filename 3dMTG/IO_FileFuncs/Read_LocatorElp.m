%Lecture des fichiers .elp format LOCATOR. Le format LOCATOR inclut le 
%format OPT3D. Si des coordonnees Brainsight sont presentes, ces coords
%seront utilisees au lieu des coordonnees LOCATOR. Toutefois, les coords
%LOCATOR seront converties de mm a m. Donc, la sortie est en metres (m).
function [matPolarisElems, matFids] = Read_LocatorElp( FullFileName )
    
    %Ouverture du fichier .elp
    fid_elp=fopen(FullFileName);
    
    %Initialisation
    matPolarisElems = {};
    matFids = [];
    
    if (fid_elp < 2)
        disp('Erreur ouverture fichier .elp');
        return;
    end
    
    %Initialisation
    matFids = zeros(3,3);
    bBrainsightLocatorFile = false;
    NbFids = 0;
    NbElems = 0;
    NbElemsColonne = 0;
    matElemsLbl = '';
    
    strline=fgets(fid_elp);

    %Boucle de lecture de lignes du fichier
    while( isempty( find( strline == -1 ) ) );

        %Gestion des commentaires Brainsight: foramt LOCATOR
        if( numel(strline) >= 12 )
            if( ~isempty( findstr( strline, sprintf( '// Polaris->' ) ) ) )
                bBrainsightLocatorFile = true;
            
            elseif( ~isempty( findstr( strline, sprintf( '//\tNZ' ) ) ) )
                stuff=strread(strline, '%s' );
                iStrtNum = strfind(strline, 'NZ')+2;
                CoordsString = strline(iStrtNum:length(strline));
                matFids(1,:) = strread( CoordsString )/1000; %mm to m
                
            elseif( ~isempty( findstr( strline, sprintf( '//\tLPA' ) ) ) )
                stuff=strread(strline, '%s' );
                iStrtNum = strfind(strline, 'LPA')+3;
                CoordsString = strline(iStrtNum:length(strline));
                matFids(2,:) = strread( CoordsString )/1000; %mm to m
                
            elseif( ~isempty( findstr( strline, sprintf( '//\tRPA' ) ) ) )
                stuff=strread(strline, '%s' );
                iStrtNum = strfind(strline, 'RPA')+3;
                CoordsString = strline(iStrtNum:length(strline));
                matFids(3,:) = strread( CoordsString )/1000; %mm to m
            end
        end
        
        if( numel(strline) >= 2 )
   
            stuff=strread(strline, '%s' );
            
            %Position of fiducials
            if( strncmp(char(stuff(1)), '%F', 2) && NbFids < 3 )
                NbFids = NbFids+1;
                CoordsString = strline(4:length(strline));
                matFids(NbFids, :) = strread( CoordsString );
                
            %Elem Type
            elseif( strncmp(char(stuff(1)), '%S', 2) )
                NbElems = NbElems+1;
                mat_tmpPolarisElems(NbElems,1) = { str2num( cell2mat(stuff(2))) };
                
            %Elem Name
            elseif( strncmp(char(stuff(1)), '%N', 2) && NbElems )
                
                strLbl = sprintf( '%-10s', char(stuff(2)) );
                strLbl = strLbl(1:10);
                matElemsLbl(NbElems,:) = strLbl;

                mat_tmpPolarisElems(NbElems,2) = { char(stuff(2)) };

                %Si c'est un fichier LOCATOR special fait par Brainsight:
                if( bBrainsightLocatorFile )
                    %Passer directement a la ligne suivante:
                    Dump = fgets(fid_elp);
                    strline = fgets(fid_elp);
                    iStrtNum = strfind(strline, char(stuff(2)))+numel(char(stuff(2)));
                    CoordsString = strline(iStrtNum:length(strline));
                    Coords=strread(CoordsString);
                    
                    %Transfert des coordonnees dans les cellules
                    for( i=1:3 ) % 3 for x, y, z
                        mat_tmpPolarisElems(NbElems,2+i) = { Coords(i)/1000 };
                    end
                %Si c'est un fichier LOCATOR normal (ex:OPT3D):
                else
                    %Coordonnées X,Y,Z,RI du trou
                    CoordsString = fgets(fid_elp);
                    Coords=strread(CoordsString);
                    
                    %Transfert des coordonnees dans les cellules
                    for( i=1:3 ) % 3 for x, y, z
                        mat_tmpPolarisElems(NbElems,2+i) = { Coords(i) };
                    end
                end
            end
        end
        
        %Lecture de la ligne suivante
        strline=fgets(fid_elp);
    end
    
    
    
    %Boucle d'effacement d'items redondands, et de moyennage des X1 X2 X3
    %X4
    NbElems = 0;
    bX1Found = false;
    bX2Found = false;
    bX3Found = false;
    bX4Found = false;
    for( iProbe = 1:size(mat_tmpPolarisElems,1) )
        if( size(matElemsLbl,2) < 10 )
            matPolarisElems = {};
            break;
        end
        %Cas d'exception: moyenner les X1
        if( ~isempty(strfind( matElemsLbl(iProbe,:), 'X1 ' )) && ~bX1Found)
            disp('X1 found');
            iOtherElems = find( matElemsLbl(:,1)  == matElemsLbl(iProbe,1) ...
                              & matElemsLbl(:,2)  == matElemsLbl(iProbe,2) ...
                              & matElemsLbl(:,3)  == matElemsLbl(iProbe,3) ...
                              & matElemsLbl(:,4)  == matElemsLbl(iProbe,4) ...
                              & matElemsLbl(:,5)  == matElemsLbl(iProbe,5) ...
                              & matElemsLbl(:,6)  == matElemsLbl(iProbe,6) ...
                              & matElemsLbl(:,7)  == matElemsLbl(iProbe,7) ...
                              & matElemsLbl(:,8)  == matElemsLbl(iProbe,8) ...
                              & matElemsLbl(:,9)  == matElemsLbl(iProbe,9) ...
                              & matElemsLbl(:,10) == matElemsLbl(iProbe,10) );
            bX1Found = true;
            NbElems = NbElems+1;
            
            matPolarisElems(NbElems,1) = mat_tmpPolarisElems(iProbe,1); %Elem Type
            matPolarisElems(NbElems,2) = mat_tmpPolarisElems(iProbe,2); %Elem Name
            matPolarisElems(NbElems,3) = { mean( [mat_tmpPolarisElems{iOtherElems,3}] ) }; %mean x
            matPolarisElems(NbElems,4) = { mean( [mat_tmpPolarisElems{iOtherElems,4}] ) }; %mean y
            matPolarisElems(NbElems,5) = { mean( [mat_tmpPolarisElems{iOtherElems,5}] ) }; %mean z
        %Cas d'exception: moyenner les X2
        elseif( ~isempty(strfind( matElemsLbl(iProbe,:), 'X2 ' )) && ~bX2Found)
            disp('X2 found');
            iOtherElems = find( matElemsLbl(:,1)  == matElemsLbl(iProbe,1) ...
                              & matElemsLbl(:,2)  == matElemsLbl(iProbe,2) ...
                              & matElemsLbl(:,3)  == matElemsLbl(iProbe,3) ...
                              & matElemsLbl(:,4)  == matElemsLbl(iProbe,4) ...
                              & matElemsLbl(:,5)  == matElemsLbl(iProbe,5) ...
                              & matElemsLbl(:,6)  == matElemsLbl(iProbe,6) ...
                              & matElemsLbl(:,7)  == matElemsLbl(iProbe,7) ...
                              & matElemsLbl(:,8)  == matElemsLbl(iProbe,8) ...
                              & matElemsLbl(:,9)  == matElemsLbl(iProbe,9) ...
                              & matElemsLbl(:,10) == matElemsLbl(iProbe,10) );
            bX2Found = true;
            NbElems = NbElems+1;
            
            matPolarisElems(NbElems,1) = mat_tmpPolarisElems(iProbe,1); %Elem Type
            matPolarisElems(NbElems,2) = mat_tmpPolarisElems(iProbe,2); %Elem Name
            matPolarisElems(NbElems,3) = { mean( [mat_tmpPolarisElems{iOtherElems,3}] ) }; %mean x
            matPolarisElems(NbElems,4) = { mean( [mat_tmpPolarisElems{iOtherElems,4}] ) }; %mean y
            matPolarisElems(NbElems,5) = { mean( [mat_tmpPolarisElems{iOtherElems,5}] ) }; %mean z
        %Cas d'exception: moyenner les X3
        elseif( ~isempty(strfind( matElemsLbl(iProbe,:), 'X3 ' )) && ~bX3Found)
            disp('X3 found');
            iOtherElems = find( matElemsLbl(:,1)  == matElemsLbl(iProbe,1) ...
                              & matElemsLbl(:,2)  == matElemsLbl(iProbe,2) ...
                              & matElemsLbl(:,3)  == matElemsLbl(iProbe,3) ...
                              & matElemsLbl(:,4)  == matElemsLbl(iProbe,4) ...
                              & matElemsLbl(:,5)  == matElemsLbl(iProbe,5) ...
                              & matElemsLbl(:,6)  == matElemsLbl(iProbe,6) ...
                              & matElemsLbl(:,7)  == matElemsLbl(iProbe,7) ...
                              & matElemsLbl(:,8)  == matElemsLbl(iProbe,8) ...
                              & matElemsLbl(:,9)  == matElemsLbl(iProbe,9) ...
                              & matElemsLbl(:,10) == matElemsLbl(iProbe,10) );
            bX3Found = true;
            NbElems = NbElems+1;
            
            matPolarisElems(NbElems,1) = mat_tmpPolarisElems(iProbe,1); %Elem Type
            matPolarisElems(NbElems,2) = mat_tmpPolarisElems(iProbe,2); %Elem Name
            matPolarisElems(NbElems,3) = { mean( [mat_tmpPolarisElems{iOtherElems,3}] ) }; %mean x
            matPolarisElems(NbElems,4) = { mean( [mat_tmpPolarisElems{iOtherElems,4}] ) }; %mean y
            matPolarisElems(NbElems,5) = { mean( [mat_tmpPolarisElems{iOtherElems,5}] ) }; %mean z%elseif( strfind( matElemsLbl(iProbe,:), 'X2 ' ) && ~bX2Found)
        
        %Cas d'exception: moyenner les X4
        elseif( ~isempty(strfind( matElemsLbl(iProbe,:), 'X4 ' )) && ~bX4Found)
            disp('X4 found');
            iOtherElems = find( matElemsLbl(:,1)  == matElemsLbl(iProbe,1) ...
                              & matElemsLbl(:,2)  == matElemsLbl(iProbe,2) ...
                              & matElemsLbl(:,3)  == matElemsLbl(iProbe,3) ...
                              & matElemsLbl(:,4)  == matElemsLbl(iProbe,4) ...
                              & matElemsLbl(:,5)  == matElemsLbl(iProbe,5) ...
                              & matElemsLbl(:,6)  == matElemsLbl(iProbe,6) ...
                              & matElemsLbl(:,7)  == matElemsLbl(iProbe,7) ...
                              & matElemsLbl(:,8)  == matElemsLbl(iProbe,8) ...
                              & matElemsLbl(:,9)  == matElemsLbl(iProbe,9) ...
                              & matElemsLbl(:,10) == matElemsLbl(iProbe,10) );
            bX4Found = true;
            NbElems = NbElems+1;
            
            matPolarisElems(NbElems,1) = mat_tmpPolarisElems(iProbe,1); %Elem Type
            matPolarisElems(NbElems,2) = mat_tmpPolarisElems(iProbe,2); %Elem Name
            matPolarisElems(NbElems,3) = { mean( [mat_tmpPolarisElems{iOtherElems,3}] ) }; %mean x
            matPolarisElems(NbElems,4) = { mean( [mat_tmpPolarisElems{iOtherElems,4}] ) }; %mean y
            matPolarisElems(NbElems,5) = { mean( [mat_tmpPolarisElems{iOtherElems,5}] ) }; %mean z%elseif( strfind( matElemsLbl(iProbe,:), 'X2 ' ) && ~bX2Found)
                                            
         %Autres cas: Prendre la derniere occurence d'un element redondant
         elseif( isempty( find(  matElemsLbl(1:iProbe-1,1)  == matElemsLbl(iProbe,1) ...
                              & matElemsLbl(1:iProbe-1,2)  == matElemsLbl(iProbe,2) ...
                              & matElemsLbl(1:iProbe-1,3)  == matElemsLbl(iProbe,3) ...
                              & matElemsLbl(1:iProbe-1,4)  == matElemsLbl(iProbe,4) ...
                              & matElemsLbl(1:iProbe-1,5)  == matElemsLbl(iProbe,5) ...
                              & matElemsLbl(1:iProbe-1,6)  == matElemsLbl(iProbe,6) ...
                              & matElemsLbl(1:iProbe-1,7)  == matElemsLbl(iProbe,7) ...
                              & matElemsLbl(1:iProbe-1,8)  == matElemsLbl(iProbe,8) ...
                              & matElemsLbl(1:iProbe-1,9)  == matElemsLbl(iProbe,9) ...
                              & matElemsLbl(1:iProbe-1,10) == matElemsLbl(iProbe,10) ) ) )
            NbElems = NbElems+1;
            
            for( iInfo=1:5 )
                matPolarisElems(NbElems,iInfo) = mat_tmpPolarisElems(iProbe,iInfo);
            end
        end
    end

    fclose(fid_elp);
end
