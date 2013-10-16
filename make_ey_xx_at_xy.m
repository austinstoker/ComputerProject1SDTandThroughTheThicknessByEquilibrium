function [ey_xx] = make_ey_xx_at_xy(Ymn,xList,yList,nList,mList,a,b)
%This will take the orientation and properties of a laminate,
% and return the ABD Matrix
ey_xx=zeros(numel(xList),numel(yList));
for i=1:numel(xList)
    x=xList(i);
    for j=1:numel(yList)
        y=yList(j);
        for i2=1:numel(mList)
            m=mList(i2);
            for j2=1:numel(nList)
                n=nList(j2);
                ey_xx(i,j)=ey_xx(i,j)+Ymn(i2,j2)*((n*pi/b)*(m*pi/a)^2)*(-1)*sin(m*pi*x/a)*sin(n*pi*y/b);
            end
        end
    end
end



