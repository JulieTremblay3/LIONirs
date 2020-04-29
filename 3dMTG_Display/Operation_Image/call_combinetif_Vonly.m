 [name,path]= uigetfile('.tif','multiselect','on')
 path = 'F:\figures_NM\'
name = {'BN_DA_Gauche.tif';
    'BN_DA_Droit.tif';
    'BN_DA_soustr.tif';}
 
 [pathout,fileout] = combinetiffileV_only(path,name)
 
  path = 'F:\figures_NM\'
name = {'BN_DF_Gauche.tif';
    'BN_DF_Droit.tif';
    'BN_DF_soustr.tif';}
 
 [pathout,fileout] = combinetiffileV_only(path,name)
 
path = 'F:\figures_NM\'
name = {'martine_DA_Gauche.tif';
    'martine_DA_Droit.tif';
    'martine_DA_soustr.tif';}
 
 [pathout,fileout] = combinetiffileV_only(path,name) 
 
 
path = 'F:\figures_NM\'
name = {'martine_DF_Gauche.tif';
    'martine_DF_Droit.tif';
    'martine_DF_soustr.tif';}
 
 [pathout,fileout] = combinetiffileV_only(path,name)