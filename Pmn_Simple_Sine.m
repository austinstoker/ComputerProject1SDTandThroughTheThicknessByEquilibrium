function [PmnSimpleSine] = Pmn_Simple_Sine(m,n,P0)
%This will take the orientation and properties of a laminate,
% and return the ABD Matrix
if m~=1 || n~=1
    PmnSimpleSine=0;
    return
else
    PmnSimpleSine=P0;
    return
end




