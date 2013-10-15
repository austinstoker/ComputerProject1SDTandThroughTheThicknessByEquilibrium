function [ T ] = T(Theta)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
Trad = (Theta*pi)/180;
m = cos(Trad);
n = sin(Trad);
T = [m^2, n^2,2*m*n; n^2,m^2,-2*n*m;-m*n,m*n, m^2-n^2];

end

