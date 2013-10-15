%This program uses Shear deformation theory to find the stresses and
%strains on a plate.  Navier solutions are used which assume Simply
%supported edges on a rectangular plate of a symmetric cross-ply laminate
%Austin Stoker Oct 2013
close all
clearvars
clc

% Set flags, counts, limits etc
doPlots=true;
PmnFunc=@Pmn_Uniform; %Which loading to use
P0=20;
numInfSum=20; %how many n,m to use
mList=1:numInfSum;
nList=1:numInfSum;
a=1; %length of the sides a
b=1; %length of sides b
xList=0:.1:1;
yList=0:.1:1;
zList=-.002:.0001:.002;


% Write the material properties
            %E1 , E2   ,  G12, v12,        ,       ,    ,v23 ,    
Mat_Types=[155E9,12.1E9,4.4E9,.248,-.018E-6,24.3E-6,.248,.458,24.3E-6];
dlmwrite('materials.txt',Mat_Types);

Angles=[90,0,0,90];
Materials=[1,1,1,1,];
Thicknesses=[.001,.001,.001,.001,];

NL=size(Angles,2);
Big(1:NL,1)=Materials(1:NL);
Big(1:NL,2)=Angles(1:NL);
Big(1:NL,3)=Thicknesses(1:NL);

dlmwrite('laminate.txt',Big);

Mat_Props(1:NL,:)=Mat_Types(Materials(1:NL),:);

% Get the expanded ABD_SDT
ABD=make_ABD_SDT(Mat_Props,Angles,Thicknesses);

% Set up mn stuff
S_hat=make_Shat(6/5,ABD,mList,a,nList,b);
b0=make_b0(S_hat);
b1=make_b1(S_hat);
b2=make_b2(S_hat);
bmn=make_bmn(S_hat,b0,b1,b2);
Pmn=make_Pmn(PmnFunc,mList,nList,P0);
Wmn=make_Wmn(Pmn,b0,bmn);
Xmn=make_Xmn(Pmn,b1,bmn);
Ymn=make_Ymn(Pmn,b2,bmn);

% Get xy stuff

P_at_xy=make_P_At_xy(Pmn,xList,yList,nList,mList,a,b);
w_at_xy=make_w_At_xy(Wmn,xList,yList,nList,mList,a,b);
Phix_at_xy=make_Phix_At_xy(Xmn,xList,yList,nList,mList,a,b);
Phiy_at_xy=make_Phiy_At_xy(Ymn,xList,yList,nList,mList,a,b);
ex_at_xy=make_ex_At_xy(Xmn,xList,yList,nList,mList,a,b);
ey_at_xy=make_ey_At_xy(Ymn,xList,yList,nList,mList,a,b);
Gamxy_at_xy=make_Gamxy_At_xy(Xmn,Ymn,xList,yList,nList,mList,a,b);
Gamxz_at_xy=make_Gamxz_At_xy(Wmn,Xmn,xList,yList,nList,mList,a,b);
Gamyz_at_xy=make_Gamyz_At_xy(Wmn,Ymn,xList,yList,nList,mList,a,b);



%% test loading
if doPlots==true
    surf(xList,yList,P_at_xy);
end

