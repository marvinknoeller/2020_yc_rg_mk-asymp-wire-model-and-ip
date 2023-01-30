function [points,t_spline,coef,br,t_stuetz] = splinepoints(p,M,t_stuetz)
if nargin ==2
    t_stuetz = zeros(1,size(p,2));
    t_stuetz(1) = 0;
    for k = 1:size(p,2)-1
       t_stuetz(k+1) = norm(p(:,k+1) - p(:,k)) + t_stuetz(k);
    end
t_stuetz = t_stuetz./t_stuetz(end);
end
y1 = p(1,:);
y2 = p(2,:);
y3 = p(3,:);
t_spline = [];
for k = 1:size(p,2)-1
    step = (t_stuetz(k+1) - t_stuetz(k))/(1*(M-1));
    t_spline = [t_spline, t_stuetz(k) : step : t_stuetz(k+1)];  
    t_spline = t_spline(1,1:end-1);     %delete the last one. So there are no double points
end

t_spline(end+1) = t_stuetz(end);
%% This is for not a knot
f1 = spline(t_stuetz,y1,t_spline);
f2 = spline(t_stuetz,y2,t_spline);
f3 = spline(t_stuetz,y3,t_spline);

pp1 = spline(t_stuetz,y1);
pp2 = spline(t_stuetz,y2);
pp3 = spline(t_stuetz,y3);
%%
coef(:,:,1) = pp1.coefs;
coef(:,:,2) = pp2.coefs;
coef(:,:,3) = pp3.coefs;
br(1,:) = pp1.breaks;
br(2,:) = pp2.breaks;
br(3,:) = pp3.breaks;
points = [f1; f2; f3];