function myconnectogram(adjacencyMatrix,label,colorMatrix,colorMap,idlist)

%PLOT NOTE 

t = linspace(-pi,pi,length(adjacencyMatrix))
factor = 1.1
coordinate = zeros(length(adjacencyMatrix),2)
for i=1:length(adjacencyMatrix)  
  plot(cos(t(i)),sin(t(i)),'o','color','k');
  edge(i)=text(cos(t(i))*factor, sin(t(i))*factor,label{i}) ; 
  if abs(t(i)) > pi/2
    edge(i).Rotation=(180*(t(i)/pi + 1));
    edge(i).HorizontalAlignment = 'right';
  else
   	edge(i).Rotation=(180*t(i)/pi+1);
   
  end
   edge(i).Color = [0 0 0];
end 
[row,col,v] = find(adjacencyMatrix);
for i=1:numel(row)
neworder(i) = colorMatrix(row(i),col(i));
end
[val,idchrono] = sort(neworder); 


for iorder = 1:length(idchrono)
     i = idchrono(iorder);
  u  = [cos(t(row(i)));sin(t(row(i)))];
            v  = [cos(t(col(i)));sin(t(col(i)))];
            x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
            y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
            r  = sqrt(x0^2 + y0^2 - 1);
            thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
            thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
% theta = [linspace(max(thetaLim),pi,50),...
%                        linspace(-pi,min(thetaLim),50)].';
    if u(1) >= 0 && v(1) >= 0 
        theta = [linspace(max(thetaLim),pi,50),...
                           linspace(-pi,min(thetaLim),50)].';
    else
       theta = linspace(thetaLim(1),thetaLim(2)).';
    end
     %DRAW ARC OF CIRCLE
    if  colorMatrix(row(i),col(i))==0
        newcolor = [0,0,0]
    else
        newcolor = colorMap(colorMatrix(row(i),col(i)),:)
    end
    line( r*cos(theta)+x0,...
              r*sin(theta)+y0,...
              'LineWidth', 4,...
              'Color', newcolor,'displayname',['id',num2str(idlist(row(i))),' ',num2str(idlist(col(i))), 
              'uicontextmenu',handles]);
end