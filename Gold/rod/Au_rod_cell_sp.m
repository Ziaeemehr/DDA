% Au-NanoRod cellular
% 007 Sloot & Drain 010 Nanosphere
clear
tic
% Simulating Sphere Nanoparticle ==========================================
[rv,N,d0]=geometrynanorod(5,20,100);    %#input:[d,L,N] 
%==========================================================================
step=15;                                %# step in wavelength
landa=(450:step:1000)*1e-9;  K=2*pi./landa; 
mm=Au_nk_sizrod;
m=mm(1,250:step:800)+i*mm(6,250:step:800);
% =========================================================================
b1=-1.891531; b2=0.1648469; b3=-1.7700004; 
% dir=[1 0 0]; pol=[0 1 0]; S=0;          %# direction:x , polarization:y
% for i=1:3
%     S=S+(dir(i)*pol(i))^2;
% end
S=1/5;
alphacm=3*d0^3*(m.^2-1)./((4*pi)*(m.^2+2));
alphaLDR=alphacm./(1+(alphacm/(d0^3)).*((b1+m.^2*b2+m.^2*b3*S).*(K*d0).^2-2/3*1i*(K*d0).^3));
a=1./alphaLDR;
% =========================================================================
V=N*d0^3; aeq=(3*V/(4*pi))^(1/3);
r=sqrt(sum(rv,2).^2);
% =========================================================================
dx = repmat(rv(:,1), 1, N) - repmat(rv(:,1).', N, 1); 
dy = repmat(rv(:,2), 1, N) - repmat(rv(:,2).', N, 1); 
dz = repmat(rv(:,3), 1, N) - repmat(rv(:,3).', N, 1);
d=cat(3,dx,dy,dz);
nd=sqrt(dx.^2+dy.^2+dz.^2);                     %# Norm of rv vectors
nx=d(:,:,1)./nd; ny=d(:,:,2)./nd; nz=d(:,:,3)./nd;
% n=cat(3,nx,ny,nz);                              %# Unit vectors
clear d
%==========================================================================
I=eye(3);
e0=1;
E0y=ones(N,1);
E0z=E0y;
Cext=zeros(1,length(landa));
Qext=zeros(1,length(landa));
Cabs=zeros(1,length(landa));
Qabs=zeros(1,length(landa));
C=cell(N);
% place=cell(N);
%==========================================================================
for s=1:length(landa)
    E0x=e0*exp(1i*K(s)*rv(:,1)); %# Light incident direction:x
    %======================================================================
    p=1;
    for i=1:N
        for j=i:N
            if      i == j
                A=a(s)*eye(3);
                A=single(A);
                C{i,j}=A;
            elseif  i < j
            A=-exp(1i*K(s)*nd(i,j))/nd(i,j)*(-K(s)^2*([nx(i,j);ny(i,j);nz(i,j)]...
                *[nx(i,j) ny(i,j) nz(i,j)]-I)+(1/nd(i,j)^2-1i*K(s)/nd(i,j))...
                *(3*[nx(i,j);ny(i,j);nz(i,j)]*[nx(i,j) ny(i,j) nz(i,j)]-I));      
            A=single(A);
            C{i,j}=A;
            C{j,i}=A;
            end
        end
    end    
%==========================================================================
    B = cell2mat(C).'; 
    C=cell(N);
    b=[E0x E0y E0z];  b=permute(b,[2 1]); b=reshape(b,3*N,1); 
    P=B\b;                              %# P=inv(B)*b
    Cext(s)=(4*pi*K(s))*imag(b'*P);
    Qext(s)=Cext(s)/(pi*aeq^2);
    Cabs(s)=(4*pi*K(s))*(imag(P.'*a(s)*conj(P))-2/3*K(s)^3*P'*P);
    Qabs(s)=Cabs(s)/(pi*aeq^2);
end
Qext=abs(Qext); Qmax=max(Qext); Qext=Qext/Qmax;
Qabs=abs(Qabs); Qmax=max(Qabs); Qabs=Qabs/Qmax;
landa=landa*1e9;
plot(landa,Qext); figure(gcf)
legend('Qext')
%==========================================================================
Qext=Qext';
Qabs=Qabs';
landa=landa';
toc