function [make_ABD_SDT] = make_ABD_SDT( Mat_Type, Angles, Thickness)
%This will take the orientation and properties of a laminate,
% and return the ABD Matrix

NL=size(Angles,2);

ABD_CLT=make_ABD(Mat_Type, Angles, Thickness);
A_CLT=ABD_CLT(1:3,1:3);
B_CLT=ABD_CLT(4:6,1:3);
D_CLT=ABD_CLT(4:6,4:6);

%Create the A matrix
A_extra=zeros(2,2);
for k=1:NL;
    Q=make_Qbar_SDT(Mat_Type(k,:),Angles(k));
    A_extra=Q(4:5,4:5)*Thickness(k)+ A_extra;
end

make_ABD_SDT=zeros(8,8);
make_ABD_SDT(1:3,1:3)=A_CLT;
make_ABD_SDT(6:8,1:3)=B_CLT;
make_ABD_SDT(1:3,6:8)=B_CLT;
make_ABD_SDT(4:5,4:5)=A_extra;
make_ABD_SDT(6:8,6:8)=D_CLT;



