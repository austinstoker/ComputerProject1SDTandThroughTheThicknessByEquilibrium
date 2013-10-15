function [ Sbar ] = Sbar( MatProp,Theta)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
E1=MatProp(1);
E2=MatProp(2);
G12=MatProp(3);
nu12=MatProp(4);
Trad = (Theta*pi)/180;
m = cos(Trad);
n = sin(Trad);

T = [m^2, n^2,2*m*n; n^2,m^2,-2*n*m;-m*n,m*n, m^2-n^2];
Tinv = T^(-1);
R = [1,0,0;0,1,0;0,0,2];
Rinv = [1,0,0;0,1,0;0,0,.5];
S = Sreduced(E1,E2,G12,nu12);

SbTemp = R*Tinv*Rinv*S*T;
Sbar = SbTemp;
end

