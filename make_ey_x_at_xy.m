function [ey_x] = make_ey_x_at_xy(Ymn,xList,yList,nList,mList,a,b)
%This will take the orientation and properties of a laminate,
% and return the ABD Matrix
ey_x=zeros(numel(xList),numel(yList));
for i=1:numel(xList)
    x=xList(i);
    for j=1:numel(yList)
        y=yList(j);
        for i2=1:numel(mList)
            m=mList(i2);
            for j2=1:numel(nList)
                n=nList(j2);
                ey_x(i,j)=ey_x(i,j)+Ymn(i2,j2)*((n*pi/b)*(m*pi/a))*cos(m*pi*x/a)*sin(n*pi*y/b);
            end
        end
    end
end



