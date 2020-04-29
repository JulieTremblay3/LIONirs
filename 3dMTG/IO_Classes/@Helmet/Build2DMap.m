%**************************************************************************
% Fonction permettant de construire une cartographie 2D du casque.
%
% La cartographie est basée sur les identificateurs de trous.
%
% Les normes d'identification du laboratoire d'IO de l'Hopital
% Ste-Justine doivent être respectées afin que cette fonction 
% soit en mesure de créer la cartographie.
%
%**************************************************************************
function oHelmet = Build2DMap( oHelmet )

    vHoles = get_vHoles(oHelmet);
    
    if( numel(vHoles) < 2 )
        return;
    end
    
    v_sPrfxs(1,:) = sprintf( '%-10s', '' ); %Columns prefixes
    v_nElems(1) = 0;
    mat_sElems(1,1,:) = sprintf( '%-10s', '' ); %Columns Elems Strings
    mat_nElems = [];
    
    for( p=1:length(vHoles) )
        if( length( vHoles(p).Label ) >= 3 )
            CurStr = sprintf( '%-10s', vHoles(p).Label );
            CurStr = CurStr(1:10);
            ac1 = int16(CurStr(1));
            ac2 = int16(CurStr(2));
            ac3 = int16(CurStr(3));
            if(    ( ac1 == int16('G') || ac1 == int16('D') ) ...
                 &&( ac2 >= int16('A') && ac2 <= int16('Z') ) )
                CurNo = [];

                if( ac3 == '(' )
                    CurPrfx = sprintf( '%-10s', CurStr(1:5) );
                    nStrBegin = find( CurStr == ')' )+1;
                    nStrBegin = nStrBegin(1);
                else
                    CurPrfx = sprintf( '%-10s', CurStr(1:2) );
                    nStrBegin = 3;
                end

                strNum = CurStr(nStrBegin:length(CurStr));
                nStrEnd = find( int16(strNum) < 48 | int16(strNum) > 57 )-1;
                if( isempty(nStrEnd) )
                    nStrEnd = length(strNum);
                end

                if( nStrEnd(1) )
                    CurNo = str2num( strNum(1:nStrEnd(1)) );
                    %disp( sprintf( 'Number found %d', CurNo ) );
                else
                    %disp( 'No number found' );
                end
                
                if( ~isempty(CurNo) )
                
                    %Recherche de prefixe semblables
                    PosC = find(   v_sPrfxs(:,1) == CurPrfx(1) & v_sPrfxs(:,2) == CurPrfx(2) ...
                                 & v_sPrfxs(:,3) == CurPrfx(3) & v_sPrfxs(:,4) == CurPrfx(4) ...
                                 & v_sPrfxs(:,5) == CurPrfx(5) );

                    %Recherche de la colonne, PosC ( Dim 1 )
                    if( isempty(PosC) )
                        if( v_sPrfxs(1,:) == sprintf( '%-10s', '' ) )
                            v_sPrfxs(1,:) = CurPrfx;
                        else
                            v_sPrfxs(size(v_sPrfxs,1)+1,:) = CurPrfx;
                        end
                        PosC = size(v_sPrfxs,1);
                        v_nElems(PosC) = 0;
                    end

                    %Recherche de la ligne, PosL ( Dim 2 )
                    CurElem = sprintf( '%-10s', vHoles(p).Label );
                    v_nElems(PosC) = v_nElems(PosC)+1;
                    PosL = v_nElems(PosC);

                    %disp( sprintf( 'PosC:%d PosL:%d CurElem:%-7s CurNo:%d', PosC, PosL, CurElem, CurNo ) );

                    mat_sElems(PosC, PosL, : ) = CurElem;
                    mat_nElems(PosC, PosL ) = CurNo;
                end
            end
        end
    end
    
    if( size(v_sPrfxs,1) == 1 && isempty( v_sPrfxs(1,:) ~= ' ' ) )
        oHelmet.Type = 'Unknown';
        return;
    end
    
    
    %Tri des lignes (Semi-sorted)
    %disp( 'Tri des lignes' );
    if( size(mat_nElems,1) && size(mat_nElems,2) )
        for( i=1:size(mat_sElems,1) )
            %disp('------------------------------------------');
            
            vSortedNumbers = sort(mat_nElems(i, : ));
            vSortedNumbers = [ vSortedNumbers( find(vSortedNumbers ~= 0) ), ...
                               vSortedNumbers( find(vSortedNumbers == 0) ) ];
            %disp( vSortedNumbers );

            Pos = [];
            for( j=1:size(mat_nElems,2) )
                v_iCorresp = find( vSortedNumbers(1,j) == mat_nElems(i,:) );
                Pos(j) = v_iCorresp(1);
            end
            %disp( Pos );
            mat_sElemsSemiSorted(i,1:length(Pos),:) = mat_sElems(i,Pos,:);
        end       
    end
    
    %Tri des colonnes: [GZ-GA,DA-DZ]
    %Tri en ordre alphabetique inverse (colonnes Gx)
    v_tmp = sortrows(v_sPrfxs( find( v_sPrfxs(:,1) == 'G' ),: ));
    for( i=1:size(v_tmp,1) )
        v_sPrfxsG(i,:) = v_tmp(size(v_tmp,1)-(i-1),:);
    end
    
    %Tri en ordre alphabetique (colonnes Dx)
    v_sPrfxsD = sortrows(v_sPrfxs( find( v_sPrfxs(:,1) == 'D' ),: ));
    
    if( ~exist( 'v_sPrfxsD' ) || ~exist( 'v_sPrfxsG' ) )
        oHelmet.Type = 'Unknown';
        return;
    else
        oHelmet.Type = 'Labo IO HSJ';
    end
    
    v_sPrfxsSorted = [v_sPrfxsG;v_sPrfxsD];
    
    Pos = [];
    %Tri des colonnes d'elements selon le tri effectue sur les prefixes
    for( i=1:size(v_sPrfxsSorted,1) )
        Pos(i) = find( v_sPrfxs(:,1) == v_sPrfxsSorted(i,1) & ...
                       v_sPrfxs(:,2) == v_sPrfxsSorted(i,2) & ...
                       v_sPrfxs(:,3) == v_sPrfxsSorted(i,3) & ...
                       v_sPrfxs(:,4) == v_sPrfxsSorted(i,4) & ...
                       v_sPrfxs(:,5) == v_sPrfxsSorted(i,5) );
    end
    
    
    mat_sElemsSorted(1:length(Pos),:,:) = mat_sElemsSemiSorted(Pos,:,:);
    
    %Les lignes de matrice sont des colonnes de casque et les colonnes de
    %matrice sont des lignes de casque.
    NbLines = size(mat_sElemsSorted,2);
    NbCols = size(mat_sElemsSorted,1);
    
    %Maintenant, construire la matrice cartographique 2D du casque:
    mat_ilHoles = zeros(NbLines,NbCols) ;
    
    %Ajout SB
    %Matrice de trous de montage
    %oHelmet.Mtg_Data.matHolesMtg = zeros(NbLines,NbCols);
    
    for( p=1:length(vHoles) )
        if( vHoles(p).Type == 400 )
            lblTmp = sprintf( '%-10s', vHoles(p).Label );
            lblTmp = lblTmp(1:10); %Tronquer a 10 car
            
            pSorted = find( mat_sElemsSorted(:,:,1) == lblTmp(1) & ...
                            mat_sElemsSorted(:,:,2) == lblTmp(2) & ...
                            mat_sElemsSorted(:,:,3) == lblTmp(3) & ...
                            mat_sElemsSorted(:,:,4) == lblTmp(4) & ...
                            mat_sElemsSorted(:,:,5) == lblTmp(5) & ...
                            mat_sElemsSorted(:,:,6) == lblTmp(6) & ...
                            mat_sElemsSorted(:,:,7) == lblTmp(7) & ...
                            mat_sElemsSorted(:,:,8) == lblTmp(8) & ...
                            mat_sElemsSorted(:,:,9) == lblTmp(9) & ...
                            mat_sElemsSorted(:,:,10) == lblTmp(10) );
                            
            if( ~isempty(pSorted) )
                %Conversion de l'indice lineaire en indice 2D (lin,col)
                %Ne pas oublier l'inversion lin-col casque-matrice
                pSorted = pSorted-1;
                Line = floor( pSorted/size(mat_sElemsSorted,1) )+1;
                Col = mod(pSorted, size(mat_sElemsSorted,1) )+1;
                %                disp( sprintf( 'Label:%s Line:%d Col:%d', vHoles(p).Label, vHoles(p).Line, vHoles(p).Col ) );
                mat_ilHoles(Line, Col) = p;
                
                %Ajout SB
                %if( p <= length(oHelmet.Mtg_Data.v_HolesMtg) )
                %    oHelmet.Mtg_Data.matHolesMtg( vHoles(p).Line, vHoles(p).Col ) = oHelmet.Mtg_Data.v_HolesMtg(p);
                %end
                
            end
        end
    end
    
    oHelmet = set_vHoles( oHelmet, vHoles );
    oHelmet.mat_ilHoles = mat_ilHoles;
