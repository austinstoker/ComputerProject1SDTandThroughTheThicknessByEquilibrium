function [ey_y] = make_ey_y_at_xy(Ymn,xList,yList,nList,mList,a,b)
%This will take the orientation and properties of a laminate,
% and return the ABD Matrix
ey_y=zeros(numel(xList),numel(yList));
for i=1:numel(xList)
    x=xList(i);
    for j=1:numel(yList)
        y=yList(j);
        for i2=1:numel(mList)
            m=mList(i2);
            for j2=1:numel(nList)
                n=nList(j2);
                ey_y(i,j)=ey_y(i,j)+Ymn(i2,j2)*((n*pi/b)^2)*sin(m*pi*x/a)*cos(n*pi*y/b);
            end
        end
    end
end



