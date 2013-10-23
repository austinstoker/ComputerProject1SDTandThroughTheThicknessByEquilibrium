%This program uses Shear deformation theory to find the stresses and
%strains on a plate.  Navier solutions are used which assume Simply
%supported edges on a rectangular plate of a symmetric cross-ply laminate
%Austin Stoker Oct 2013
close all
clearvars
clc

% Set flags, counts, limits etc
mycase=1;
doPlots=true;
if mycase==0  %verification case a
    PmnFunc=@Pmn_Simple_Sine;
    P0=10000;
    numInfSum=1;
    Angles=[0,90,90,0];
    numLayers=size(Angles,2);
    a=.1;
    b=.1;
    xList=[a/2];
    yList=[b/2];
    zsteps=1000;
    xPtQuarter = a/2;
    yPtQuarter = b/2;
    t=.000127;
    
elseif mycase==1  %verification case b
    PmnFunc=@Pmn_Uniform;
    P0=10000;
    numInfSum=500;
    Angles=[0,90,90,0];
    numLayers=size(Angles,2);
    a=.1;
    b=.1;
    xList=[a/4];
    yList=[b/4];
    zsteps=1000;
    xPtQuarter = a/4;
    yPtQuarter = b/4;
    t=.000127;
    
elseif mycase==2 %a from the project assignment
    PmnFunc=@Pmn_Uniform;
    P0=10000;
    numInfSum=500;
    Angles=[0,90,90,0];
    numLayers=size(Angles,2);
    a=.1;
    b=.1;
    xList=[a/4,a/2];
    yList=[b/4,b/2];
    zsteps=1000;
    xPtQuarter = a/4;
    yPtQuarter = b/4;
    xPtCenter = a/2;
    yPtCenter = b/2;
    
    t=.000127; 
    
elseif mycase==3 %b from the project assignment
    PmnFunc=@Pmn_Uniform;
    P0=10000;
    numInfSum=500;
    Angles=[0,90,90,0,90,90,0];%
    numLayers=size(Angles,2);
    a=.2; %
    b=.1;
    xList=[a/4,a/2];
    yList=[b/4,b/2];
    zsteps=1000;
    xPtQuarter = a/4;
    yPtQuarter = b/4;
    xPtCenter = a/2;
    yPtCenter = b/2;
    t=.000127; 

elseif mycase==4 %a from the project assignment
    PmnFunc=@Pmn_Simple_Sine; %
    P0=10000;
    numInfSum=1; %
    Angles=[0,90,90,0];
    numLayers=size(Angles,2);
    a=.1;
    b=.1;
    xList=[a/4,a/2];
    yList=[b/4,b/2];
    zsteps=1000;
    xPtQuarter = a/4;
    yPtQuarter = b/4;
    xPtCenter = a/2;
    yPtCenter = b/2;
    
    t=.000127; 
    
elseif mycase==5 %d from the project assignment
    PmnFunc=@Pmn_HydrostaticVaryingInX;
    P0=10000;
    numInfSum=500;
    Angles=[0,90,90,0,90,90,0];%
    numLayers=size(Angles,2);
    a=.2; %
    b=.1;
    xList=[a/4,a/2];
    yList=[b/4,b/2];
    zsteps=1000;
    xPtQuarter = a/4;
    yPtQuarter = b/4;
    xPtCenter = a/2;
    yPtCenter = b/2;
    t=.000127; 
end
   
zList=-t*numLayers/2+t/zsteps:t/zsteps:t*numLayers/2-t/zsteps;
mList=1:numInfSum;
nList=1:numInfSum;

K=5/6;


% Write the material properties
            %E1 , E2   ,  G12, v12,        ,       ,    ,v23 ,    
Mat_Types=[155E9,12.1E9,4.4E9,.248,-.018E-6,24.3E-6,.248,.458,24.3E-6];
dlmwrite('materials.txt',Mat_Types);


Materials=ones(1,numel(Angles)); % [1,1,1,1,1,1,1,1];
Thicknesses= t*ones(1,numel(Angles)); %[t,t,t,t,t,t,t,t];

NL=size(Angles,2);
Big(1:NL,1)=Materials(1:NL);
Big(1:NL,2)=Angles(1:NL);
Big(1:NL,3)=Thicknesses(1:NL);

dlmwrite('laminate.txt',Big);

Mat_Props(1:NL,:)=Mat_Types(Materials(1:NL),:);

% Get the expanded ABD_SDT
ABD=make_ABD_SDT(Mat_Props,Angles,Thicknesses);

% Set up mn stuff
S_hat=make_Shat(K,ABD,mList,a,nList,b);
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


% Setup strains
ex_xyz=zeros(numel(xList),numel(yList),numel(zList));
ey_xyz=zeros(numel(xList),numel(yList),numel(zList));
Gamxy_xyz=zeros(numel(xList),numel(yList),numel(zList));
Strains_xyz=zeros(numel(xList),numel(yList),numel(zList),5);
for k=1:numel(zList)
    ex_xyz(:,:,k)=ex_at_xy.*-zList(k); % 1) only these three vary with z
    ey_xyz(:,:,k)=ey_at_xy.*-zList(k); % 2
    Gamxy_xyz(:,:,k)=Gamxy_at_xy.*zList(k); % 3
    %Build the Strains matrix
    Strains_xyz(:,:,k,1)=ex_xyz(:,:,k);
    Strains_xyz(:,:,k,2)=ey_xyz(:,:,k);
    Strains_xyz(:,:,k,3)=Gamxy_xyz(:,:,k);
    Strains_xyz(:,:,k,4)=Gamyz_at_xy(:,:);
    Strains_xyz(:,:,k,5)=Gamxz_at_xy(:,:);
end


% Setup Qbar for each layer
Qbar=zeros(5,5,NL);
for k=1:NL
    Qbar(:,:,k)=make_Qbar_SDT(Mat_Props(k,:),Angles(k));
end

% Get the Stresses
Stress_xyz=zeros(numel(xList),numel(yList),numel(zList),5);
for i=1:numel(xList)
    x=xList(i);
    for j=1:numel(yList)
        y=yList(j);
        for k=1:numel(zList)
            z=zList(k);
            layer=get_layer(z,Thicknesses);
            thisStrain=zeros(5,1);
            for r=1:5
                thisStrain(r)=Strains_xyz(i,j,k,r);
            end
            Stress_xyz(i,j,k,:)=Qbar(:,:,layer)*thisStrain;
        end
    end
end

% Now the fun stuff, transverse junk
% tau_xz
ex_x=make_ex_x_at_xy(Xmn,xList,yList,nList,mList,a,b);
ey_x=make_ey_x_at_xy(Ymn,xList,yList,nList,mList,a,b);
Gamxy_y=make_Gamxy_y_at_xy(Xmn,Ymn,xList,yList,nList,mList,a,b);

tau_xz_e_xyz=zeros(numel(xList),numel(yList),numel(zList));
sig_x_x = tau_xz_e_xyz;
tau_xy_y = tau_xz_e_xyz;
for i=1:numel(xList)
    for j=1:numel(yList)
        tau_xz_e_xyz(i,j,1)=0;
        for k=2:numel(zList)
            z=zList(k);
            layer=get_layer(z,Thicknesses);
            sig_x_x(i,j,k) = Qbar(1,1,layer)*ex_x(i,j)+Qbar(1,2,layer)*ey_x(i,j);
            tau_xy_y(i,j,k) = Qbar(3,3,layer)*Gamxy_y(i,j);
            tau_xz_e_xyz(i,j,k)=tau_xz_e_xyz(i,j,k-1)+(sig_x_x(i,j,k)+tau_xy_y(i,j,k))*(zList(k)^2-zList(k-1)^2)/2;
        end
    end
end

% tau_yz
ey_y=make_ey_y_at_xy(Ymn,xList,yList,nList,mList,a,b);
ex_y=make_ex_y_at_xy(Xmn,xList,yList,nList,mList,a,b);
Gamxy_x=make_Gamxy_x_at_xy(Xmn,Ymn,xList,yList,nList,mList,a,b);

tau_yz_e_xyz=zeros(numel(xList),numel(yList),numel(zList));
sig_y_y = tau_yz_e_xyz;
tau_xy_x = tau_yz_e_xyz;
for i=1:numel(xList)
    for j=1:numel(yList)
        tau_yz_e_xyz(i,j,1)=0;
        for k=2:numel(zList)
            z=zList(k);
            layer=get_layer(z,Thicknesses);
            sig_y_y(i,j,k) = Qbar(2,2,layer)*ey_y(i,j)+Qbar(1,2,layer)*ex_y(i,j);
            tau_xy_x(i,j,k) = Qbar(3,3,layer)*Gamxy_x(i,j);
            tau_yz_e_xyz(i,j,k)=tau_yz_e_xyz(i,j,k-1)+(sig_y_y(i,j,k)+tau_xy_x(i,j,k))*(zList(k)^2-zList(k-1)^2)/2;
        end
    end
end

%sigma z
ex_xx=make_ex_xx_at_xy(Xmn,xList,yList,nList,mList,a,b);
ey_xx=make_ey_xx_at_xy(Ymn,xList,yList,nList,mList,a,b);
Gamxy_xy=make_Gamxy_xy_at_xy(Xmn,Ymn,xList,yList,nList,mList,a,b);
ey_yy=make_ey_yy_at_xy(Ymn,xList,yList,nList,mList,a,b);
ex_yy=make_ex_yy_at_xy(Xmn,xList,yList,nList,mList,a,b);

T_hat=zeros(numel(xList),numel(yList),numel(zList));
for i=1:numel(xList)
    for j=1:numel(yList)
        for k=1:numel(zList)
            z=zList(k);
            layer = get_layer(z,Thicknesses);
            T_hat(i,j,k) = Qbar(1,1,layer)*ex_xx(i,j)+ Qbar(1,2,layer)*(ey_xx(i,j))... %sigma_x_xx
                          -2*(Qbar(3,3,layer)*Gamxy_xy(i,j))... %2*tauxy_xy
                          +Qbar(1,2,layer)*(ex_yy(i,j))+Qbar(2,2,layer)*ey_yy(i,j); %sigma_y_yy
        end
    end
end

sig_z_e_xyz=zeros(numel(xList),numel(yList),numel(zList));
for i=1:numel(xList)
    for j=1:numel(yList)
        sig_z_e_xyz(i,j,1)=-P0;
        for k=2:numel(zList)
            z=zList(k);
            sigz_hat_0 = 0;
            for j2=2:k-1
                sigz_hat_0 = sigz_hat_0+T_hat(j2)*(zList(j2)^2-zList(j2-1)^2)/2;
            end
            sig_z_e_xyz(i,j,k)=sig_z_e_xyz(i,j,k-1)-...
             sigz_hat_0*(zList(k)-zList(k-1))-T_hat(k)*...
             ((zList(k)^3)/6 +(zList(k-1)^3)/3 - (zList(k)*zList(k-1)^2)/2);
        end
    end
end

%%
if doPlots==true 
    [~,i] = min(abs(xList-xPtQuarter));
    [~,j] = min(abs(yList-yPtQuarter));
    
    sigx=zeros(numel(zList),1);
    sigy=zeros(numel(zList),1);
    epsx=zeros(numel(zList),1);
    tauxy=zeros(numel(zList,1));
    tauxz=zeros(numel(zList,1));
    tauyz=zeros(numel(zList,1));
    sigz=zeros(numel(zList,1));
    for k=1:numel(zList)
        sigx(k)=Stress_xyz(i,j,k,1);
        sigy(k)=Stress_xyz(i,j,k,2);
        tauxy(k)=Stress_xyz(i,j,k,3);
        epsx(k)=Strains_xyz(i,j,k,1);
        tauxz(k)=tau_xz_e_xyz(i,j,k);
        tauyz(k)=tau_yz_e_xyz(i,j,k);
        sigz(k)=sig_z_e_xyz(i,j,k);
    end
    figure(1)
    plot(sigx,zList)
    title('sig x');
    
    figure(2)
    plot(sigy,zList)
    title('sig y');
    
    figure(3)
    plot(tauxy,zList)
    title('tau_x_y')
    
    figure(4)
    plot(tauxz,zList);
    title('tau_x_z')
    
    figure(5)
    plot(tauyz,zList);
    title('tau_y_z')
    
    figure(6)
    plot(sigz,zList);
    title('sig_z')
    
    figure(7)
    %surf(xList,yList,w_at_xy);
    
end



%% test loading
if doPlots==true
    figure(7);
    %surf(xList,yList,P_at_xy);
end

