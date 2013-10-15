function [ Qbar_SDT ] = make_Qbar_SDT(MatProp, Theta )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

m=cosd(Theta);
n=sind(Theta);

E2=MatProp(2);
v23=MatProp(8);
G12=MatProp(3);
G13=G12;

G23=E2/(2*(1+v23));
C44=G23;
C55=G13;

Q44=m^2*C44+n^2*C55;
Q45=m*n*(C55-C44);
Q55=n^2*C44+m^2*C55;
Qbar_extra=[Q44, Q45; Q45,Q55];

Qbar_SDT(1:3,1:3)=make_Qbar(MatProp,Theta);
Qbar_SDT(4:5,4:5)=Qbar_extra;



end

