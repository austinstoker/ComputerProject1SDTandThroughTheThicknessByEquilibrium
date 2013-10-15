function [Gamyz_xy] = make_Gamyz_At_xy(Wmn,Ymn,xList,yList,nList,mList,a,b)
%This will take the orientation and properties of a laminate,
% and return the ABD Matrix
Gamyz_xy=zeros(numel(xList),numel(yList));
for i=1:numel(xList)
    x=xList(i);
    for j=1:numel(yList)
        y=yList(j);
        for i2=1:numel(mList)
            m=mList(i2);
            for j2=1:numel(nList)
                n=nList(j2);
                Gamyz_xy(i,j)=Gamyz_xy(i,j)+Wmn(i2,j2)*n*pi/b*sin(m*pi*x/a)*cos(n*pi*y/b)...
                                           +Ymn(i2,j2)*sin(m*pi*x/a)*cos(n*pi*y/b);
            end
        end
    end
end



