% Etablir l'intervalle d'affichage des identificateurs de trous
%
% Intrants: oDispOpt  : Objet IO_DisplayOptions
%           Interval  : Entier, Intervalle d'affichage des identificateurs
%
% Extrants: oDispOpt : Objet modifie
%
%**************************************************************************
function oDispOpt = set_HoleIdInterval( oDispOpt, Interval )

    oDispOpt.Lbl_HelmetHoleId_Interval = Interval;
    