function [wxy] = make_w_At_xy(Wmn,xList,yList,nList,mList,a,b)
%This will take the orientation and properties of a laminate,
% and return the ABD Matrix
wxy=zeros(numel(xList),numel(yList));
for i=1:numel(xList)
    x=xList(i);
    for j=1:numel(yList)
        y=yList(j);
        for i2=numel(mList)
            m=mList(i2);
            for j2=1:numel(nlist)
                n=nList(j2);
                wxy(i,j)=wxy(i,j)+Wmn(i2,j2)*sin(m*pi*x/a)*sin(n*pi*y/b);
            end
        end
    end
end



