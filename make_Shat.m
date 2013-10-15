function [ make_Shat ] = make_Shat( K,ABD_SDT,mList,a,nList,b)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
A55=ABD_SDT(5,5);
A44=ABD_SDT(4,4);
D11=ABD_SDT(6,6);
D12=ABD_SDT(6,7);
D22=ABD_SDT(7,7);
D66=ABD_SDT(8,8);

make_Shat=zeros(numel(mList),numel(nList));
for i=1:numel(mList)
    m=mList(i);
    for j=1:numel(nList)
        n=nList(j);
        alm=m*pi/a;
        betn=n*pi/b;
        S11=K*(A55*alm^2+A44*betn^2);
        S12=K*A55*alm;
        S13=K*A44*betn;
        S22=K*A55+D11*alm^2+D66*betn^2;
        S23=(D12+D66)*alm*betn;
        S33=K*A44+D22*betn^2+D66*alm^2;
        make_Shat(i,j,1:3,1:3)=[S11,S12,S13; S12,S22,S23; S13,S23,S33];
    end
end

end

