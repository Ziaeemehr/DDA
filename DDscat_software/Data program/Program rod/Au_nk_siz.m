function result=Au_nk_siz(lambdamin,lambdamax,r,source)
% r Radius of NPs e.g: r=[3 5 9 15 21 30 40 50] in (nm)
% Wavelemgth interval  [lambdamin lambdamax]  in (nm)
% Sizing Dielectric Constants 
% Modified Refractive index of different size of gold nanoparticles
% Introducing Variables ==================================================
%wp             Plasma frequency [s^-1]
%gama_bulk      Damping constant in Drude model [s^-1]
%lambda          Wavelength of incident light (nm)
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
% Recalling Constants
wp=1.3e16;           
gama_bulk=1.64e14;   
c=299792458e9;                  %# [nm/s] 
lambda=(lambdamin:lambdamax);   %# [nm]
w=2*pi*c./lambda;    
% D=1;
vf=14.1e14;
%# Recalling Experimental bulk data
var=epsilonAu(lambdamin,lambdamax,source);   
eps=var(4,:);                   %# eps1+ i eps2
% Size Correction for Gold nanorods =======================================
for j=1:length(r)
    eps_siz=eps+wp^2./(w.*(w+i./gama_bulk))-wp^2./(w.*(w+i./gama_bulk+i*vf./(r(j))));
    eps1=real(eps_siz);
    eps2=imag(eps_siz);
    
%   Refractive index - corrected by size of the particle ==================
    n_siz(j,:)=(0.5*((eps1.^2+eps2.^2).^0.5+eps1)).^0.5;
    k_siz(j,:)=(0.5*((eps1.^2+eps2.^2).^0.5-eps1)).^0.5;
end
% Size Correction for Gold nanospheres ====================================
% for j=1:length(r)
%     gama=gama_bulk+D*vf./(r(j));
%     eps_boun=eps-1+(wp^2./(w.^2+1i*w.*gama_bulk));
%     eps_1=real(eps_boun);
%     eps_2=imag(eps_boun);
% 
%     % Dielectric function - corrected by size of the particle ===========
%     eps1_siz=1-(wp^2./(w.^2+gama.^2))+eps_1;
%     eps2_siz=(wp^2*gama)./((w.^2+gama.^2).*w)+eps_2;
% 
%     % Refractive index - corrected by size of the particle ==============
%     n_siz(j,:)=(0.5*((eps1_siz.^2+eps2_siz.^2).^0.5+eps1_siz)).^0.5;
%     k_siz(j,:)=(0.5*((eps1_siz.^2+eps2_siz.^2).^0.5-eps1_siz)).^0.5;
% end
%n_bulk=var(2,:);
%k_bulk=var(3,:);
result=[n_siz;k_siz];