function [Gamxy_xy] = make_Gamxy_At_xy(Xmn,Ymn,xList,yList,nList,mList,a,b)
%This will take the orientation and properties of a laminate,
% and return the ABD Matrix
Gamxy_xy=zeros(numel(xList),numel(yList));
for i=1:numel(xList)
    x=xList(i);
    for j=1:numel(yList)
        y=yList(j);
        for i2=1:numel(mList)
            m=mList(i2);
            for j2=1:numel(nList)
                n=nList(j2);
                Gamxy_xy(i,j)=Gamxy_xy(i,j)+Xmn(i2,j2)*n*pi/b*cos(m*pi*x/a)*cos(n*pi*y/b)...
                                           +Ymn(i2,j2)*m*pi/a*cos(m*pi*x/a)*cos(n*pi*y/b);
            end
        end
    end
end



