function result=Au_nk_sizrod()   
% for Rod Leff=4*V/S
% Sizing Dielectric Constants 
% Modified Refractive index of different size of gold nanoparticles
%% Introducing Variables ==================================================
%wp             Plasma frequency [s^-1]
%gama_bulk      Damping constant in Drude model [s^-1]
%landa          Wavelength of incident light (nm)
%w              Angular frequency of incident light [s^-1]
%c              Speed of light [nm]
%eps_boun       Dielectric function of bound-electron
%eps            Experimental values of Dielectric function of Bulk marerial
%gama           Modified Damping factor in Drude model
%D              Is a constant that includes details of the scattering
%               processes
%vf             Electron velocity at the Fermi surface [nm/s]
%r              Radius of the particle [nm]
%eps1_siz       Real part of dielectric function for a particle of radius r
%eps2_siz       Imag part of dielectric function for a particle of radius r
%n,k            Interpolated values of Experimental Complex refractive index
% =========================================================================

clear; clc;
%% Recalling Constants
r=3.75;       %# Radius of NP
L=[25 27 40 50];
%# Wavelenght 
landa=(201:1200);               %# [nm]              
n_siz=zeros(length(L),length(landa));
k_siz=n_siz;
wp=1.3e16;           
gama_bulk=1.64e14;   
c=299792458e9;                  %# [nm/s] 
w=2*pi*c./landa;    
% D=1.1;
vf=14.1e14;
var=InterpAu;                   %# Recalling Experimental bulk data
eps=var(3,:);                   %# eps1+ i eps2
Leff=4*pi*r.^2*L./(2*pi*r.*L+2*pi*r.^2);
% =========================================================================
for j=1:length(L)
    % gama=gama_bulk+D*vf./(r(j));     %# For Sphere
	gama=gama_bulk+vf./Leff(j);        %# For Rod
    eps_boun=eps-1+(wp^2./(w.^2+1i*w.*gama_bulk));
    eps_1=real(eps_boun);
    eps_2=imag(eps_boun);

    % Dielectric function - corrected by size of the particle =============
    eps1_siz=1-(wp^2./(w.^2+gama.^2))+eps_1;
    eps2_siz=(wp^2*gama)./((w.^2+gama.^2).*w)+eps_2;

    % Refractive index - corrected by size of the particle ================
    n_siz(j,:)=(0.5*((eps1_siz.^2+eps2_siz.^2).^0.5+eps1_siz)).^0.5;
    k_siz(j,:)=(0.5*((eps1_siz.^2+eps2_siz.^2).^0.5-eps1_siz)).^0.5;
end

var2=Johnson_exp;
landa_exp=var2(12,:);     
n_exp=var2(4,:);
k_exp=var2(5,:);
%% Interpolation Experimental data ========================================
n=interp1(landa_exp,n_exp,landa,'spline');    
k=interp1(landa_exp,k_exp,landa,'spline');    
%eps1_bulk=n.^2-k.^2;
%eps2_bulk=2*n.*k;
result=[n_siz;n;k_siz;k];
%result is a 18x1000 matrix 
