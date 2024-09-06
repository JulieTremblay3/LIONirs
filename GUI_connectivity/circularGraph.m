classdef circularGraph < handle
% CIRCULARGRAPH Plot an interactive circular graph to illustrate connections in a network.
%
%% Syntax
% circularGraph(X)
% circularGraph(X,'PropertyName',propertyvalue,...)
% h = circularGraph(...)
%
%% Description
% A 'circular graph' is a visualization of a network of nodes and their
% connections. The nodes are laid out along a circle, and the connections
% are drawn within the circle. Click on a node to make the connections that
% emanate from it more visible or less visible. Click on the 'Show All'
% button to make all nodes and their connections visible. Click on the
% 'Hide All' button to make all nodes and their connections less visible.
%
% Required input arguments.
% X : A symmetric matrix of numeric or logical values.
%
% Optional properties.
% Colormap : A N by 3 matrix of [r g b] triples, where N is the 
%            length(adjacenyMatrix).
% Label    : A cell array of N strings.
%%
% Copyright 2016 The MathWorks, Inc.
  properties
    Node = node(0,0); % Array of nodes
    ColorMap;         % Colormap
    Label;            % Cell array of strings
    ShowButton;       % Turn all nodes on
    HideButton;       % Turn all nodes off
    colorlink;        % link color
  end
  
  methods
    function this = circularGraph(adjacencyMatrix,varargin)
      % Constructor
      p = inputParser;
 
      defaultColorMap = parula(length(adjacencyMatrix));
      defaultLabel = cell(length(adjacencyMatrix));
      for i = 1:length(defaultLabel)
        defaultLabel{i} = num2str(i);
      end
      colmanual = ones(length(defaultLabel));
      defaultcolorlink = ones(length(defaultLabel))
      addRequired(p,'adjacencyMatrix',@(x)(isnumeric(x) || islogical(x)));
      addParameter(p,'ColorMap',defaultColorMap,@(colormap)length(colormap) == length(adjacencyMatrix));
      addParameter(p,'Label'   ,defaultLabel   ,@iscell);
      addParameter(p,'colorlink'   ,colmanual  ,@isnumeric);
      parse(p,adjacencyMatrix,varargin{:});
     
      this.ColorMap = p.Results.ColorMap;
      this.Label    = p.Results.Label;
      this.colorlink    = p.Results.colorlink
      this.ShowButton = uicontrol(...
        'Style','pushbutton',...
        'Position',[0 40 80 40],...
        'String','Show All',...
        'Callback',@circularGraph.showNodes,...
        'UserData',this,...
        'visible','off');
      
      this.HideButton = uicontrol(...
        'Style','pushbutton',...
        'Position',[0 0 80 40],...
        'String','Hide All',...
        'Callback',@circularGraph.hideNodes,...
        'UserData',this,...
        'visible','off');
      
      fig = gcf;
      set(fig,...
        'UserData',this,...
        'CloseRequestFcn',@circularGraph.CloseRequestFcn);
      
      % Draw the nodes
      delete(this.Node);
      t = linspace(-pi,pi,length(adjacencyMatrix) + 1).'; % theta for each node
      extent = zeros(length(adjacencyMatrix),1);
      for i = 1:length(adjacencyMatrix)
        this.Node(i) = node(cos(t(i)),sin(t(i)));
        this.Node(i).Color = this.ColorMap(i,:);
        this.Node(i).Label = this.Label{i};
       % set(gco,'fontsize',14)
      end 
       tmp = get(gca,'children');
       for i=1:numel(tmp)
      if strcmp(get(tmp(i),'type'),'text')
       set(tmp(i),'fontsize',14) ;
      end
       end
       
      % Find non-zero values of s and their indices
      [row,col,v] = find(adjacencyMatrix);
      
      % Calculate line widths based on values of s (stored in v).
      minLineWidth  = 3;
      lineWidthCoef = 5;
      lineWidth = v./max(v);
      if sum(lineWidth) == numel(lineWidth) % all lines are the same width.
        lineWidth = repmat(minLineWidth,numel(lineWidth),1);
      else % lines of variable width.
        lineWidth = lineWidthCoef*lineWidth + minLineWidth;
      end
      
      % Draw connections on the Poincare hyperbolic disk.
      %
      % Equation of the circles on the disk:
      % x^2 + y^2 
      % + 2*(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1))*x 
      % - 2*(u(1)-v(1))/(u(1)*v(2)-u(2)*v(1))*y + 1 = 0,
      % where u and v are points on the boundary.
      %
      % Standard form of equation of a circle
      % (x - x0)^2 + (y - y0)^2 = r^2
      %
      % Therefore we can identify
      % x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
      % y0 = (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
      % r^2 = x0^2 + y0^2 - 1
      colorid =  jet(100); %JT
      %colorid = flipud(colorid);
      for i=1:numel(row)
      neworder(i) = this.colorlink(row(i),col(i));
      end
      [val,idchrono] = sort(neworder);      
      
  
      for iorder = 1:length(v)
          i = idchrono(iorder);
        if row(i) ~= col(i)
          newcolor =  colorid(this.colorlink(row(i),col(i)),:) ;

          if abs(row(i) - col(i)) - length(adjacencyMatrix)/2 == 0 
            % points are diametric, so draw a straight line
            u = [cos(t(row(i)));sin(t(row(i)))];
            v = [cos(t(col(i)));sin(t(col(i)))];
            
            %colormanual
            %newcolor = this.ColorMap(row(i),:);
            this.Node(row(i)).Connection(end+1) = line(...
              [u(1);v(1)],...
              [u(2);v(2)],...
              'LineWidth', lineWidth(i),...
              'Color', newcolor,...
              'PickableParts','none');
          else % points are not diametric, so draw an arc
            u  = [cos(t(row(i)));sin(t(row(i)))];
            v  = [cos(t(col(i)));sin(t(col(i)))];
            x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
            y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
            r  = sqrt(x0^2 + y0^2 - 1);
            thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
            thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
            
            if u(1) >= 0 && v(1) >= 0 
              % ensure the arc is within the unit disk
              theta = [linspace(max(thetaLim),pi,50),...
                       linspace(-pi,min(thetaLim),50)].';
            else
              theta = linspace(thetaLim(1),thetaLim(2)).';
            end
           
            this.Node(row(i)).Connection(end+1) = line(...
              r*cos(theta)+x0,...
              r*sin(theta)+y0,...
              'LineWidth', lineWidth(i),...
              'Color', newcolor,...
              'PickableParts','none');
          end
        end
      end
      
      axis image;
      ax = gca;
      for i = 1:length(adjacencyMatrix)
        extent(i) = this.Node(i).Extent;
      end
      extent = max(extent(:));
      ax.XLim = ax.XLim + extent*[-1 1];
      fudgeFactor = 1.75; % Not sure why this is necessary. Eyeballed it.
      ax.YLim = ax.YLim + fudgeFactor*extent*[-1 1];
      ax.Visible = 'off';
      ax.SortMethod = 'depth';
      
      fig = gcf;
      fig.Color = [1 1 1];
    end
    
  end
  
  methods (Static = true)
    function showNodes(this,~)
      % Callback for 'Show All' button
      n = this.UserData.Node;
      for i = 1:length(n)
        n(i).Visible = true;
      end
    end
    
    function hideNodes(this,~)
      % Callback for 'Hide All' button
      n = this.UserData.Node;
      for i = 1:length(n)
        n(i).Visible = false;
      end
    end
    
    function CloseRequestFcn(this,~)
      % Callback for figure CloseRequestFcn
      c = this.UserData;
      for i = 1:length(c.Node)
        delete(c.Node(i));
      end
      delete(gcf);
    end
    
  end
  
end