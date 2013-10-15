function [ Sreduced ] = Sreduced( E1,E2,G12,nu12)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Stemp = [1/E1,-nu12/E1,0;-nu12/E1 ,1/E2,0 ;0,0,1/G12];
Sreduced=Stemp;

end

