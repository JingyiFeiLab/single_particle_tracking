function [X,Y] = calcEllipse(varargin) 
% function [X,Y] = calculateEllipse(x, y, a, b, angle, steps) 
%# This functions returns points to draw an ellipse 
%# 
%# @param x X coordinate 
%# @param y Y coordinate 
%# @param a Semimajor axis 
%# @param b Semiminor axis 
%# @param angle Angle of the ellipse (in rad) 
%# 
% Source: http://stackoverflow.com/questions/2153768/draw-ellipse-and-ellipsoid-in-matlab/24531259#24531259 
% Modified by Christian Fässler

steps = 360; 

if nargin == 1 || nargin == 2 
x = varargin{1}.X0_in; 
y = varargin{1}.Y0_in; 
a = varargin{1}.a; 
b = varargin{1}.b;
angle = varargin{1}.angleToX;
%angle = varargin{1}.phi; 
if nargin == 2 
steps = varargin{2}; 
end 
else if nargin == 5 || nargin == 6 
x = varargin{1}; 
y = varargin{2}; 
a = varargin{3}; 
b = varargin{4}; 
angle = varargin{5}; 
if nargin == 6 
steps = varargin{6}; 
end 
else 
error('Wrong input'); 
end 
end 

beta = -angle*(pi/180);
%beta = -angle; 
sinbeta = sin(beta); 
cosbeta = cos(beta);

alpha = linspace(0,2*pi, 2*steps)'; 
% sinalpha = (sin(alpha)).^(.7); 
% cosalpha = (cos(alpha)).^(.7);

sinalpha = sin(alpha); 
cosalpha = cos(alpha);

if a>=b
    a = a+1;
    X = sign(cosalpha).*a.*abs(cosalpha).^(.5) * cosbeta - sign(sinalpha).*b.*abs(sinalpha).^(.5)*sinbeta;
    Y = sign(sinalpha).*b.*abs(sinalpha).^(.5)*cosbeta + sign(cosalpha).*a.*abs(cosalpha).^(.5) * sinbeta;
    X = X+x;
    Y = Y+y;
%     X = abs(round(x + (a * sign(cosalpha).*abs(cosalpha).^(.666) * cosbeta - b * sign(sinalpha).*abs(sinalpha).^(.666) * sinbeta))); 
%     Y = abs(round(y + (a * sign(cosalpha).*abs(cosalpha).^(.666) * sinbeta + b * sign(sinalpha).*abs(sinalpha).^(.666) * cosbeta)));
%     X = round(x + (a * cosalpha * cosbeta - b * sinalpha * sinbeta)); 
%     Y = round(y + (a * cosalpha * sinbeta + b * sinalpha * cosbeta));
    X = int32(X);
    Y = int32(Y);
elseif a <b
    
    X = sign(cosalpha).*b.*abs(cosalpha).^(.5) * cosbeta - sign(sinalpha).*a.*abs(sinalpha).^(.5)*sinbeta;
    Y = sign(sinalpha).*a.*abs(sinalpha).^(.5)*cosbeta + sign(cosalpha).*b.*abs(cosalpha).^(.5) * sinbeta;
    X = X+x;
    Y = Y+y;
%     X = round(x + (b * cosalpha * cosbeta - a * sinalpha * sinbeta)); 
%     Y = round(y + (b * cosalpha * sinbeta + a * sinalpha * cosbeta));
    X = int32(X);
    Y = int32(Y);
end

for i = 1:length(X)-1
    if abs(X(i+1)-X(i)) > 1 || abs(Y(i+1)-Y(i)) > 1
        line = makeLine([X(i+1),X(i)],[Y(i+1),Y(i)]);
        X = [X; line(:,2)];
        Y = [Y; line(:,1)];
    end
end
    



if nargout==1
    X = [int32(X) int32(Y)]; 
end 


end
