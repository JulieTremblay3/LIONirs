% Permet de savoir l'intervalle d'affichage des identificateurs de trous
%
% Intrants: oDispOpt  : Objet IO_DisplayOptions
%
% Extrants: Interval : Entier, Intervalle d'affichage des identificateurs
%
%**************************************************************************
function Interval = get_HoleIdInterval( oDispOpt )

    Interval = oDispOpt.Lbl_HelmetHoleId_Interval;