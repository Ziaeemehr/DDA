function [rv,Numberdipole,d0]=geometrySphere(R,N)
% ÂÏÑÓ Ïåí ÏæŞØÈí åÇ ÈÑÇí í˜ ˜Ñå Èå ÔÚÇÚ ãÚáæã ˜å ÈÇ N Ïæ ŞØÈí ÊŞÑíÈ ÒÏå
% ÔÏå ÇÓÊ.                                                            
% R                # Radius of nanoparticle in meter
% N                # Number of dipoles 
% d0               # Distance between dipoles
%                  # Input is approximat but output is exact
% result           # Coordinates of dipoles
clc
%# help geometrySphere
R=R*1e-9;
d0=(4*pi/(3*N))^(1/3)*R;
% d0=d0*1e-9;        %# in nanometer
% R=R*1e-9;          %# in nanometer
x=-R:d0:R; 
% n=floor(2*R/d0);
% x=linspace(-R,R,n);
y=x; z=x;
% Main Loop ===============================================================
t=1;
for i=1:length(x)
	for j=1:length(y)
		for k=1:length(z)
			u=x(i)^2+y(j)^2+z(k)^2;
			if u < R^2
                rv(t,:)=[x(i) y(j) z(k)];
                t=t+1;
			end
		end
	end
end

Numberdipole=length(rv);
fprintf('N=%4.0f d0= %5.1e \n', Numberdipole, d0)

% Ploting the results =====================================================
%  plot3(rv(:,1),rv(:,2),rv(:,3),'o')
% xlabel('x'); ylabel('y'); zlabel('z');
