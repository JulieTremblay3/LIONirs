function ChannelLabels = ConvertmlIDsrs2label(NIRS)
%Look at identification source number to convert in a source label
%according to the device information 
NIRS.Cf.dev.n = 'NIRx'; %HARD CODE FOR ELAN
    for id=1:size(NIRS.Cf.H.C.id,2) 
         if strcmp(NIRS.Cf.dev.n,'NIRx') %NIRx
            srs = SDPairs2strboxy(NIRS.Cf.H.C.id(2,id));
            det = SDDet2strboxy(NIRS.Cf.H.C.id(3,id));
            ChannelLabels{id,1} = [srs, '_', det];
         else                             %On ajoute la nomenclature Ste-Justine 
            srs = SDPairs2strboxy_ISS(NIRS.Cf.H.C.id(2,id));
            det = SDDet2strboxy_ISS(NIRS.Cf.H.C.id(3,id));
            ChannelLabels{id,1} = [srs, '_', det];
        end
    end