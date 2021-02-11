function result=epsilonAu(lambdamin,lambdamax,source)
% Data range: [0.125 -- 11.2] (eV) or [9918--110] (nm)
%DEFINITION DES CONSTANTES'
% hb=1.05457266E-34;    %# h/(2pi) J.s
c=299792458e9;      
h_e=4.1356992e-15;      %# eV.s

if          strcmp(source,'palik')
%PALIK
load Au_Palik           
    a=Au_Palik;
    E_exp = a(:,1);
    n_exp = a(:,2);
    k_exp = a(:,3);
elseif      strcmp(source,'JC')
%J-C
load Au_JC           
    a=Au_JC;
    E_exp = a(:,1);
    n_exp = a(:,2);
    k_exp = a(:,3);
end
lambda_exp=h_e*c./E_exp;
lambda=(lambdamin:lambdamax);
%INTERPOLATION OF THE DIELECTRIC FUNCTION
n=(interp1(lambda_exp,n_exp,lambda,'cubic'));
k=(interp1(lambda_exp,k_exp,lambda,'cubic'));
eps1=n.^2-k.^2;
eps2=2*n.*k;
epsilon=eps1+1i*eps2;
result=[lambda;n;k;epsilon];