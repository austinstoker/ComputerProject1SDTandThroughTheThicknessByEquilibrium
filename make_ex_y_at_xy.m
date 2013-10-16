function [ex_y] = make_ex_y_at_xy(Xmn,xList,yList,nList,mList,a,b)
%This will take the orientation and properties of a laminate,
% and return the ABD Matrix
ex_y=zeros(numel(xList),numel(yList));
for i=1:numel(xList)
    x=xList(i);
    for j=1:numel(yList)
        y=yList(j);
        for i2=1:numel(mList)
            m=mList(i2);
            for j2=1:numel(nList)
                n=nList(j2);
                ex_y(i,j)=ex_y(i,j)+Xmn(i2,j2)*(m*pi/a)*(n*pi/b)*sin(m*pi*x/a)*cos(n*pi*y/b);
            end
        end
    end
end



