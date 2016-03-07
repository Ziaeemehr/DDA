clear
clc
[rv,N,d0]=geometrySphere(5,2000);
tic
dx1 = meshgrid(rv(:,1)) - meshgrid(rv(:,1))' ;
dy1 = meshgrid(rv(:,2)) - meshgrid(rv(:,2))' ;
dz1 = meshgrid(rv(:,3)) - meshgrid(rv(:,3))' ;
dx1=dx1.'; dy1=dy1.'; dz1=dz1.';
toc

% tic
% for i=1:N
%     for j=1:N
%         dx(i,j)=rv(i,1)-rv(j,1);
%         dy(i,j)=rv(i,2)-rv(j,2);
%         dz(i,j)=rv(i,3)-rv(j,3);
%     end
% end
% toc
tic 
dx2 = repmat(rv(:,1), 1, N) - repmat(rv(:,1).', N, 1); 
dy2 = repmat(rv(:,2), 1, N) - repmat(rv(:,2).', N, 1); 
dz2 = repmat(rv(:,3), 1, N) - repmat(rv(:,3).', N, 1); 
toc