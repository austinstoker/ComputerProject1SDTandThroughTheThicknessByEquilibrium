function [ Qbar ] = make_Qbar( MatProp,Theta )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Qbar=make_Sbar(MatProp,Theta)^-1;


end

