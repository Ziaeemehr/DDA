function result=InterpAu()
% Interpolation of refractive index for gold nanoparticles with precision
% 1(nm) by Spline interpolation

% Programmer            Date
%==========================================================================
% A.Ziaee Mehr          07.89
%==========================================================================

%% Introducing Variables

%var1           Experimental Valus of Johnson and Christy 
%E              Energy of incident light (eV)
%h_e            h/e(in eV)
%c              Speed of light (in nm)
%landa_exp      Wavelength of incident light on surface
%n_exp          Real part of refractive index of Au
%k_exp          Imag part of refractive index of Au
%landa          Wavelength for interpolation points
%n              Interpolated values for n
%k              Interpolated values for k
%epsilon        Dielectric function = epsilon1+i epsilon2
%Re_epsilon     %epsilon1
%Img_epsilon    %epsilon2

%% Initialisation =========================================================
clear
clc
%% E in eV experimental ===================================================
var1=Johnson_exp;
landa_exp=var1(12,:);

%% Gold optical constants Experimental ====================================
n_exp=var1(4,:);
k_exp=var1(5,:);

%% Interpolation Spline
landa=201:1200;         %# Wavelength [nm] 
n=interp1(landa_exp,n_exp,landa,'spline');    
k=interp1(landa_exp,k_exp,landa,'spline');    
% Dielectric function  epsilon= epsilon1+i epsilon2
Re_epsilon=n.^2-k.^2;               
Img_epsilon=2*n.*k;      
epsilon=Re_epsilon+1i*Img_epsilon;   
result=[Re_epsilon;Img_epsilon;epsilon;n;k];
%% ========================================================================
% result is a 5x1000 matrix 
% Re_epsilon        1
% Img_epsilon       2
% epsilon           3
% n                 4
% k                 5
% For a fast check
% plot(landa_exp,n_exp,'o',landa,n,landa_exp,k_exp,'*',landa,k)