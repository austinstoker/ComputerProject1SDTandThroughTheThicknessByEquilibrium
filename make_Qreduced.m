function [ Qreduced ] = make_Qreduced( E1,E2,G12,nu12)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Qreduced=Sreduced( E1,E2,G12,nu12 )^-1;

end

