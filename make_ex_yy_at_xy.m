function [ex_yy] = make_ex_yy_at_xy(Xmn,xList,yList,nList,mList,a,b)
%This will take the orientation and properties of a laminate,
% and return the ABD Matrix
ex_yy=zeros(numel(xList),numel(yList));
for i=1:numel(xList)
    x=xList(i);
    for j=1:numel(yList)
        y=yList(j);
        for i2=1:numel(mList)
            m=mList(i2);
            for j2=1:numel(nList)
                n=nList(j2);
                ex_yy(i,j)=ex_yy(i,j)+Xmn(i2,j2)*((m*pi/a)*(n*pi/b)^2)*sin(m*pi*x/a)*(-1)*sin(n*pi*y/b);
            end
        end
    end
end



