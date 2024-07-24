% Samyukta Ramnath
% BITS Pilani K.K. Birla Goa Campus
% March 2017
%function: cone. Generates a right circular cone with axis orientation
%aperture angle and height specified.
% inputs : theta: aperture angle of cone
%          d    : vector describing orientation of axis of cone. Format:
%          2x3 vector with x,y,z coordinates of two points lying on the
%          axis of the cone. d(1,:) = (x1,y1,z1); d(2,:) = (x2,y2,z2).
%          h    : height of cone       
function [X3,Y3,Z3,X,Y,Z] = cone(theta,d,h)
%    theta = 45;
    
    r = h*tan(pi*theta/180);
    m = h/r;
    [R,A] = meshgrid(linspace(0,r,11),linspace(0,2*pi,41));
    % Generate cone about Z axis with given aperture angle and height
    X = R .* cos(A);
    Y = R .* sin(A);
    Z = m*R;
    % Cone around the z-axis, point at the origin
    % find coefficients of the axis vector xi + yj + zk
    x = d(2,1)-d(1,1);
    y = d(2,2)-d(1,2);
    z = d(2,3)-d(1,3);
    
    % find angle made by axis vector with X axis
    phix = atan2(y,x);
    % find angle made by axis vector with Z axis
    phiz = atan2(sqrt(x^2 + y^2),(z));
    
    % Rotate once about Z axis 
    X1 = X*cos(phiz)+Z*sin(phiz);
    Y1 = Y;
    Z1 = -X*sin(phiz)+Z*cos(phiz);
    % Rotate about X axis
    X3 = X1*cos(phix)-Y1*sin(phix);
    Y3 = X1*sin(phix)+Y1*cos(phix);
    Z3 = Z1;
    
end