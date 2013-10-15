function [Pmn_HydrostaticVaryingInX] = Pmn_HydrostaticVaryingInX(m,n,P0)
%This will take the orientation and properties of a laminate,
% and return the ABD Matrix
if  mod(n,2) == 0 %if n is even
    Pmn_HydrostaticVaryingInX=0;
    return
else %if n is odd
    Pmn_HydrostaticVaryingInX=8*P0/(pi^2*m*n)*(-1)^(m+1);
    return
end




