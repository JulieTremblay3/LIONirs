function handles=plot_axes_linkMAP(varargin)
%plot_axes_linkConnectogramm(handles, newfigure) 

%plot_axes_linkConnectogramm(handles, newfigure) 
handles = varargin{1};
if numel(varargin)>1; newfigure = varargin{2};else; newfigure = 0;end 
if newfigure ==0    
    axes(handles.axes_viewlink);
    cla;
    hold on;
else
    figure;hold on;
end
 
 load(get(handles.edit_linkSettingmap2dmap,'string'))
 image(map(:,:,1:3))