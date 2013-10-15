function [Pmn_Uniform] = Pmn_Uniform(m,n,P0)
%This will take the orientation and properties of a laminate,
% and return the ABD Matrix
if mod(n,2) == 0 || mod(m,2) == 0 %if n or m is even
    Pmn_Uniform=0;
    return
else %if m and n are odd
    Pmn_Uniform=16*P0/(pi^2*m*n);
    return
end




