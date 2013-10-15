function [ make_Pmn ] = make_Pmn(PmnFunc,mList,nList,P0)
%UNTITLED3 Summary of this function goes here
%   Detaild explanation goes here
make_Pmn=zeros(numel(mList),numel(nList));

for i=1:numel(mList)
    m=mList(i);
    for j=1:numel(nList)
        n=nList(j);
        make_Pmn(i,j)=PmnFunc(m,n,P0);
    end
end
end

