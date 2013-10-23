function [ABD] = make_ABD( Mat_Type, Angles, Thickness)
%This will take the orientation and properties of a laminate,
% and return the ABD Matrix



NL=size(Angles,2); %number of layers



%Total thickness of the laminate
Tt=0;
for k=1:NL
    Tt=Thickness(k) + Tt;
end


%positions of lamina edges, z(0) in the book = z(1) here because
%matlab doesn't do 0 index
z(1)=0-Tt/2;
for k=1:NL
   z(k+1)=Thickness(k) + z(k);
end


%Create the A matris
A=zeros(3,3);
for k=1:NL;
    A=make_Qbar(Mat_Type(k,:),Angles(k))*Thickness(k)+ A;
end

MinA=max(max(A))/10000;
for j=1:3
    for k=1:3
        if abs(A(j,k))<=MinA
            A(j,k)=0;
        end
    end
end


%The B Matrix
B=zeros(3,3);
for k=1:NL
    B=(.5)*make_Qbar(Mat_Type(k,:),Angles(k))*(z(k+1)^2 - z(k)^2) + B;
end
%zero out the tiny values of B
for k=1:3
    for i=1:3
        if abs(10E16*B(i,k))< abs(max(max(A)))
            B(i,k)=0;
        end
    end
end


%The D Matrix
D=zeros(3,3);
for k=1:NL
    D=(1/3)*make_Qbar(Mat_Type(k,:),Angles(k))*(z(k+1)^3 - z(k)^3) + D;
end
MinD=max(max(D))/10000;
for j=1:3
    for k=1:3
        if abs(D(j,k))<=MinD
            D(j,k)=0;
        end
    end
end


%Populate the ABD
ABD=zeros(6,6);
for i=1:3
    for j=1:3
        ABD(i,j)=A(i,j);
    end
    for j=4:6
        ABD(i,j)=B(i,j-3);
    end
end 
for i=4:6
    for j=1:3
        ABD(i,j)=B(i-3,j);
    end
    for j=4:6
        ABD(i,j)=D(i-3,j-3);
    end
end 
