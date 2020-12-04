function [X,Y,Z,nx,ny,nz] = actin_polarity_fit_bundle_plane(x,y,z,Imdim)
%%%%%%%%%%
%%%%%%%%%%
% This program fits a plane to filament center of mass coordinates.
%
% INPUT
% x --- x coordinates of the sampling points /in pixel
% y --- y coordinates of the sampling points /in pixel
% z --- z coordinates of the sampling points /in pixel
% Imdim --- Image dimension of the sampling points /in pixel
%
% OUTPUT
% X,Y,Z --- x,y,z coordinates of points element of fitted plane

% Plane fit
[A,B,C] = plane_fit(x,y,z);

% Plane calculation
[Xn,Yn] = meshgrid(1:1:Imdim);
Y = single(Xn);
X = single(Yn);
Z = (A.*X)+(B.*Y)+C;

% Hessian normal form
nx = -A./sqrt((-A).^2+(-B).^2+(-C).^2);
ny = -B./sqrt((-A).^2+(-B).^2+(-C).^2);
nz = +1./sqrt((-A).^2+(-B).^2+(-C).^2);




function [A,B,C]=plane_fit(x,y,z)

% function [A,B,C]=plane_fit(x,y,z)
% ------------------------------------------------------------------------
%   Fit a plane to x,y,z data.
%   [A,B,C]=plane_fit(x,y,z) calculates the coefficients A,B,C that fit the data
%   defined by the vectors x,y,z. 
%
%   Uses command svd
%
%   %EXAMPLE: 
%
%         [x,y]=meshgrid(linspace(0,10,20),linspace(0,10,20));
%         a=1; b=2; c=-2;
%         z=(a*x)+(b*y)+c;
%         x=x(:); y=y(:); z=z(:);
%         z=z+(randn(length(z),1));
%         [A,B,C]=plane_fit(x,y,z); 
%         [X,Y]=meshgrid(linspace(min(x),max(x),20),linspace(min(y),max(y),20));
%         Z=(A*X)+(B*Y)+C;
%         plot3(x,y,z,'r.'); hold on; grid on;
%         surf(X,Y,Z,'FaceColor','g'); alpha(0.5);
%         title(['a=',num2str(a), ', A=',num2str(A),', b=',num2str(b),', B=',num2str(B),', c=',num2str(c),', C=',num2str(C)]);
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 14/08/2008
% ------------------------------------------------------------------------

P=[mean(x),mean(y),mean(z)];
[U,S,V]=svd([x-P(1),y-P(2),z-P(3)],0);
N=-1/V(end,end)*V(:,end);
A=N(1); B=N(2); C=-P*N;

