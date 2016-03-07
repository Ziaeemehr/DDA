clear , clc
Ax=[1 2;3 4];
Ay=[5 6;7 8];
Az=[9 10;11 12];
% A=cat(3,Ax,Ay,Az)
% B=permute(A,[3 2 1])
N=2;
% Ax=rand(N); Ay=rand(N); Az=rand(N);
A = cat(3, Ax, Ay, Az);   % [N-by-N-by-3]
B=permute(A,[3 2 1]);
B=reshape(B,3,N^2);
F = zeros(3, 3, N^2);
for i = 1:3,  
	for j = 1:3,    
		Ai = A(:,:,i); 
		Aj = A(:,:,j);  
		P = Ai * Aj; 
		F(i,j,:) = reshape(P, [1, 1, N^2]); 
	end
 end
 %# then we can write
%  B = mat2cell(F,3,3,ones(N^2,1));
%  B = reshape(B,N,N)'; 
%  B = cell2mat(B);
 