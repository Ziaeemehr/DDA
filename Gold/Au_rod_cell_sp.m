% Au-NanoRod cellular
% 007 Sloot & Drain 010 Nanosphere
clear
tic
% Simulating Sphere Nanoparticle ==========================================
[rv,N,d0]=geometrynanorod(7.5,25,270);    %#input:[d,L,N] 
dir=[1 0 0]; pol=[0 1 0];               %# direction , polarization
%==========================================================================
step=10;                                %# step in wavelength
landa=(450:step:1200)*1e-9;  K=2*pi./landa; 
%mm=Au_nk_sizrod;
mm=Au_nk_sizrod_L_fix;
m=mm(3,250:step:1000)+i*mm(8,250:step:1000);
% =========================================================================
b1=-1.891531; b2=0.1648469; b3=-1.7700004; S=0;
for i=1:3
    S=S+(dir(i)*pol(i))^2;
end
% S=1/5;
alphacm=3*d0^3*(m.^2-1)./((4*pi)*(m.^2+2));
alphaLDR=alphacm./(1+(alphacm/(d0^3)).*((b1+m.^2*b2+m.^2*b3*S).*(K*d0).^2-2/3*1i*(K*d0).^3));
a=1./alphaLDR; a=single(a);
% =========================================================================
V=N*d0^3; aeq=(3*V/(4*pi))^(1/3);
r=sqrt(sum(rv,2).^2);
dx = repmat(rv(:,1), 1, N) - repmat(rv(:,1).', N, 1); 
dy = repmat(rv(:,2), 1, N) - repmat(rv(:,2).', N, 1); 
dz = repmat(rv(:,3), 1, N) - repmat(rv(:,3).', N, 1);
nd=sqrt(dx.^2+dy.^2+dz.^2);   %# Norm of rv vectors
nx=dx./nd; ny=dy./nd; nz=dz./nd;
%==========================================================================
I=eye(3);
e0=1;
E0y=ones(N,1);
E0z=E0y;
Cext=zeros(1,length(landa));
Qext=zeros(1,length(landa));
Cabs=zeros(1,length(landa));
Qabs=zeros(1,length(landa));
diam=find(eye(N));              % find index of main diagonal
for s=1:length(landa)
    E0x=e0*exp(1i*K(s)*rv(:,1)); %# Light incident direction:x
    C = cell(N);
    for  i=1:N
        for j=i:N
            if i==j
                continue
            else
                n=[nx(i,j) ny(i,j) nz(i,j)]; ndij=nd(i,j);
                nn=n(:)*n; 
                A=-exp(1i*K(s)*nd(i,j))/ndij*(-K(s)^2*(nn-I)+...
                    (1/ndij^2-1i*K(s)/ndij)*(3*nn-I)); 
                A=single(A);
                C{i,j} = A; C{j,i} = A;
            end
        end
    end
    C(diam)={a(s)*eye(3);};
    C = cell2mat(C);  
    b=[E0x E0y E0z];  b=permute(b,[2 1]); b=reshape(b,3*N,1); 
    P=C\b;                              %# P=inv(B)*b
    Cext(s)=(4*pi*K(s))*imag(b'*P);
    Qext(s)=Cext(s)/(pi*aeq^2);
    Cabs(s)=-(4*pi*K(s))*(imag(P.'*a(s)*conj(P))-2/3*K(s)^3*(P'*P));
    Qabs(s)=Cabs(s)/(pi*aeq^2);
end
% Qext=abs(Qext); Qmax=max(Qext); Qext=Qext/Qmax;
% Qabs=abs(Qabs); Qmax=max(Qabs); Qabs=Qabs/Qmax;
landa=landa*1e9;
plot(landa,Qext,landa,Qabs); figure(gcf)
legend('Qext','Qabs')
%==========================================================================
Qext=Qext';
Qabs=Qabs';
landa=landa';
toc