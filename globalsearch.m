function [x] = globalsearch(ax,by)
%GLOBALSEARCH 确定壁面点位置
%   此处显示详细说明
YYY= evalin('base', 'YYY');
xmin= evalin('base', 'xmin');
xmax= evalin('base', 'xmax');
MIN=@(x) (x-ax).^2+(ppval(YYY,x)-by).^2;
rng default
opts = optimoptions(@fmincon,'Algorithm','sqp','MaxIterations',500);
problem=createOptimProblem('fmincon','objective',MIN,'x0',ax,'lb',ax-0.1,'ub',ax+0.1,'options',opts);
gs=GlobalSearch;
[x,~]=run(gs,problem);
end
