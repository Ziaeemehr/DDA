function [rv,Numberdipole,d0]=geometrynanorod(d,L,N)
% d                # diameter of Nanorod
% N                # approximate Number of dipoles 
%                  # Input is approximat but output is exact
% d0               # lattice spacing
% result           # Coordinates of dipoles
% =========================================================================
d=d*1e-9; L=L*1e-9; 
R=d/2;
d0=(4*pi/(3*N))^(1/3)*R;
x=-R:d0:R; 
y=-L/2:d0:L/2; z=x;
% Main Loop ===============================================================
R2=R^2;
t=1;
for i=1:length(x)
	for j=1:length(y)
		for k=1:length(z)
			u=x(i)^2+z(k)^2;
			if u < R2
                rv(t,:)=[x(i) y(j) z(k)];
                t=t+1;
			end
		end
	end
end

Numberdipole=length(rv);
fprintf('N=%4.0f d0= %5.1e \n', Numberdipole, d0)

% Ploting the results =====================================================
% plot3(rv(:,1),rv(:,2),rv(:,3),'o')
% xlabel('x'); ylabel('y'); zlabel('z');  
