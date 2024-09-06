function obj = set_Helmet(obj, oHelmet)
% IO_Project_Data/get_Helmet 
% Returns the Helmet object
if( isa(oHelmet, 'Helmet') )
    obj.m_Helmet = oHelmet;
else
    warndlg( 'Arg not a Helmet' );
    oHelmet
end
