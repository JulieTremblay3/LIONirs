function nirs_average_defaults
%
% Sets the defaults for NIRS

global nirsHSJ

nirsHSJ.dataprocessing.averaging.savenirs = false;
nirsHSJ.dataprocessing.averaging.triggervalue = 1;
nirsHSJ.dataprocessing.averaging.pretime = 5;
nirsHSJ.dataprocessing.averaging.posttime = 30;
nirsHSJ.dataprocessing.averaging.averagetype = 1;
nirsHSJ.dataprocessing.averaging.badintervalratio = 0.6;