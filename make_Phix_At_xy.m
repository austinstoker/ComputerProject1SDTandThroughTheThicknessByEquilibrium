function [Phix_xy] = make_Phix_At_xy(Xmn,xList,yList,nList,mList,a,b)
%This will take the orientation and properties of a laminate,
% and return the ABD Matrix
Phix_xy=zeros(numel(xList),numel(yList));
for i=1:numel(xList)
    x=xList(i);
    for j=1:numel(yList)
        y=yList(j);
        for i2=1:numel(mList)
            m=mList(i2);
            for j2=1:numel(nList)
                n=nList(j2);
                Phix_xy(i,j)=Phix_xy(i,j)+Xmn(i2,j2)*cos(m*pi*x/a)*sin(n*pi*y/b);
            end
        end
    end
end



