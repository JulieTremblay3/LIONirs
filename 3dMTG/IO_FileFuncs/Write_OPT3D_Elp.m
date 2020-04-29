%Creation d'un fichier Elp format compatible OPT3D.
%Les positions des sources et des detecteurs sont ecrites (le montage).
%
%NOTE: Toutes les positions sont prolongees le long de leur normales 
%      d'un longueur equivalente l'epaisseur de la peau. Si aucun MRI ou
%      aucune segmentation n'est attachee au projet, alors la profondeur
%      est zero par defaut, ou peut etre une epaisseur constante selon le
%      choix defini dans l'ecran de configuration de parametres du montage.
function Write_OPT3D_Elp( oHelmet, PathFileName,all,fitonskin,flag_MRIfid )

    if( ~PathFileName ) 
        return;
    end

    v_Holes = get_vHoles(oHelmet);
    sMtg = get_Mtg(oHelmet);

    if isfield(sMtg,'v_pEle')
        v_pSrcDet = [sMtg.v_pDet, sMtg.v_pSrc, sMtg.v_pEle];
    else
        v_pSrcDet = [sMtg.v_pDet, sMtg.v_pSrc];
    end
    oHelmet = calc_Center( oHelmet);
    center = get_Center( oHelmet )
    %Replace le centre à 0,0,0     
    fid  = sMtg.matFiducials;
%     if size(fid,1)==6
%         fid = fid' 
%     end
%     if fid(1,2)==0 & fid(1,3)==0 & fid(2,3)==0 & fid(3,3)== 0
%         rotate_ok = 0;
%     else 
%         rotate_ok = 1;
%     end
%     rotate_ok=0;
%     if rotate_ok
%         %centre Rpa Lpa
%         fid(:,4) = [1;1;1];
%         cx = [fid(2,1)+( fid(3,1)-fid(2,1))/2, fid(2,2)+( fid(3,2)-fid(2,2))/2 ,fid(2,3)+( fid(3,3)-fid(2,3))/2 ]
%         cm = makehgtform('translate',[-cx(1) -cx(2) -cx(3)])
%         fidc = fid*cm';
%         alphanas = atan(fidc(1,3)/fidc(1,2))
%         mRpa = makehgtform('xrotate', -alphanas)
%         fidc = fidc * mRpa';       
%         theta = atan(fidc(1,2)/fidc(1,1));
%         mRpt = makehgtform('zrotate',-theta);
%         fidc = fidc * mRpt';      
%         thetaax = atan(fidc(2,3)/fidc(2,2));
%         mLpa = makehgtform('xrotate',-thetaax);
%         fidc = fidc * mLpa';      
%         sMtg.matFiducials = fidc(1:3,1:3);  
%     end
%   makehgtform('xrotate',  
%     figure
%     hold on
%     coord = [];
%     for i = 1:numel(v_Holes)
%         plot3(v_Holes(i).Coord.x,v_Holes(i).Coord.y,v_Holes(i).Coord.z,'*')
%         coord = [coord; v_Holes(i).Coord.x,v_Holes(i).Coord.y,v_Holes(i).Coord.z,1];
%     end
%     
    %Il y a eu un déplacement manuel de la référence entre la tête et les
    %trous
%     if size(sMtg.matFiducials,1)==6
%         offset_all = (sMtg.matFiducials(1:3,:)-sMtg.matFiducials(4:6,:))
%         offset_dx = max(offset_all(:,1));
%         offset_dy = max(offset_all(:,2));
%         offset_dz = max(offset_all(:,3)); %NON UTILISÉ
%     else
%         offset_dx = 0;
%         offset_dy = 0;
%         offset_dz = 0;
%     end

%     if 0
%         fid(:,4) = [1;1;1];
%         %cx = [fid(2,1)+( fid(3,1)-fid(2,1))/2, fid(2,2)+(
%         %fid(3,2)-fid(2,2))/2 ,fid(2,3)+( fid(3,3)-fid(2,3))/2 ]
%         cx = [0,0,0,0];
%         cm = makehgtform('translate',[-cx(1) -cx(2) -cx(3)])
%         fidc = fid*cm';
% %         alphanas = atan(fidc(1,3)/fidc(1,2))
%         alphanas = 0;
%         mRpa = makehgtform('xrotate', -alphanas)
%         fidc = fidc * mRpa';              
% %         theta = atan(fidc(1,2)/fidc(1,1));
%         theta = 0;
%         mRpt = makehgtform('zrotate',-theta);
%         fidc = fidc * mRpt';      
%         thetaax  = pi/24; %rotation lpa rpa
% %         thetaax = atan(fidc(2,3)/fidc(2,2));
%         mLpa = makehgtform('xrotate',-thetaax);
%         fidc = fidc * mLpa';      
%         sMtg.matFiducials = fidc(1:3,1:3);
%         rotate_ok = 1;
%     end
%     

    fid_elp = fopen(PathFileName, 'w');
    if fid_elp == -1
        return
    end
        
    fprintf(fid_elp, '3\t2\n');
    fprintf(fid_elp, '//Probe file\n');
    fprintf(fid_elp, '//Minor revision number\n');
    fprintf(fid_elp, '1\n');
    fprintf(fid_elp, '//ProbeName\n');
    fprintf(fid_elp, '%%%%N\tName\n');
    fprintf(fid_elp, '//Probe type, number of sensors\n');
    if all
       fprintf(fid_elp, '0\t%d\n', size(v_Holes,2))
    else
        fprintf(fid_elp, '0\t%d\n', numel(v_pSrcDet));
    end
    fprintf(fid_elp, '//Position of fiducials X+, Y+, Y- on the subject\n');
%     if flag_MRIfid
%         fprintf(fid_elp, '%%F\t%.6f\t%.6f\t%.6f\n', sMtg.matFiducials(4,1), sMtg.matFiducials(4,2), sMtg.matFiducials(4,3) );
%         fprintf(fid_elp, '%%F\t%.6f\t%.6f\t%.6f\n', sMtg.matFiducials(5,1), sMtg.matFiducials(5,2), sMtg.matFiducials(5,3) );
%         fprintf(fid_elp, '%%F\t%.6f\t%.6f\t%.6f\n', sMtg.matFiducials(6,1), sMtg.matFiducials(6,2), sMtg.matFiducials(6,3) );
%     else
    	fprintf(fid_elp, '%%F\t%.6f\t%.6f\t%.6f\n', sMtg.matFiducials(1,1), sMtg.matFiducials(1,2), sMtg.matFiducials(1,3) );
        fprintf(fid_elp, '%%F\t%.6f\t%.6f\t%.6f\n', sMtg.matFiducials(2,1), sMtg.matFiducials(2,2), sMtg.matFiducials(2,3) );
        fprintf(fid_elp, '%%F\t%.6f\t%.6f\t%.6f\n', sMtg.matFiducials(3,1), sMtg.matFiducials(3,2), sMtg.matFiducials(3,3) );
%     end
%     
    Counter = 0;
    if all
         for p = 1 : size(v_Holes,2)
            Counter = Counter+1;
%             if rotate_ok
%                 c(1) = v_Holes(p).Coord.x;
%                 c(2) = v_Holes(p).Coord.y;
%                 c(3) = v_Holes(p).Coord.z;
%                 c(4) = 1;
%                 ct = c*cm';
%                 ct = ct*mRpa';
%                 ct = ct*mRpt';
%                 ct = ct* mLpa';
%                 v_Holes(p).Coord.x =ct(1);
%                 v_Holes(p).Coord.y =ct(2);
%                 v_Holes(p).Coord.z = ct(3);
%             end
            fprintf(fid_elp, '//Sensor type\n');
            fprintf(fid_elp, '%%S\t%d\n', v_Holes(p).Type );
            fprintf(fid_elp, '//Sensor name and data for sensor# %d\n', Counter );
            fprintf(fid_elp, '%%N\t%s\n', v_Holes(p).Label );
            if fitonskin
                fprintf(fid_elp, '%.6f\t%.6f\t%.6f\n', v_Holes(p).Coord.x-v_Holes(p).Normal.x*v_Holes(p).SkinDepth, ...
                                                   v_Holes(p).Coord.y-v_Holes(p).Normal.y*v_Holes(p).SkinDepth, ...
                                                   v_Holes(p).Coord.z-v_Holes(p).Normal.z*v_Holes(p).SkinDepth );
            else
                fprintf(fid_elp, '%.6f\t%.6f\t%.6f\n', v_Holes(p).Coord.x, ...
                                                       v_Holes(p).Coord.y, ...
                                                       v_Holes(p).Coord.z );
            end
        end
    else
           for( Pos=1:numel(v_pSrcDet) )
            p = v_pSrcDet(Pos);
            Counter = Counter+1;
            if rotate_ok
                c(1) = v_Holes(p).Coord.x;
                c(2) = v_Holes(p).Coord.y;
                c(3) = v_Holes(p).Coord.z;
                c(4) = 1;
                ct = c*cm';
                ct = ct*mRpa';
                ct = ct*mRpt';
                ct = ct* mLpa';            
                v_Holes(p).Coord.x =ct(1);
                v_Holes(p).Coord.y =ct(2);
                v_Holes(p).Coord.z = ct(3);
            end
            fprintf(fid_elp, '//Sensor type\n');
            fprintf(fid_elp, '%%S\t%d\n', v_Holes(p).Type );
            fprintf(fid_elp, '//Sensor name and data for sensor# %d\n', Counter );
            fprintf(fid_elp, '%%N\t%s\n', v_Holes(p).Label );
                if fitonskin
                     fprintf(fid_elp, '%.6f\t%.6f\t%.6f\n', v_Holes(p).Coord.x-v_Holes(p).Normal.x*v_Holes(p).SkinDepth, ...
                                                       v_Holes(p).Coord.y-v_Holes(p).Normal.y*v_Holes(p).SkinDepth, ...
                                                       v_Holes(p).Coord.z-v_Holes(p).Normal.z*v_Holes(p).SkinDepth );
                else
                    fprintf(fid_elp, '%.6f\t%.6f\t%.6f\n', v_Holes(p).Coord.x, ...
                                                       v_Holes(p).Coord.y, ...
                                                       v_Holes(p).Coord.z );
                end
           end
    end
    fclose(fid_elp);
    
   