clear
N=300;
Ax=rand(N); Ay=rand(N); Az=rand(N);     %# [N-by-N]
% tic
% t=1;
% F=zeros(3,3,N^2);
% for  i=1:N
%     for j=1:N
%         F(:,:,t)= [Ax(i,j);Ay(i,j);Az(i,j)]*[Ax(i,j) Ay(i,j) Az(i,j)];
%         t=t+1;                          %# t is just a counter
%     end
% end
% %# then we can write
% B = mat2cell(F,3,3,ones(N^2,1));
% B = reshape(B,N,N)'; 
% B = cell2mat(B);
% toc
% =========================================================================
tic
B = cell(N);
for index = 1:N^2
  A = [Ax(index) Ay(index) Az(index)];
  B{index} = A(:)*A;
end
B = cell2mat(B);
toc