close all
clear all
clc


myCase=1; %1 Uniform, 2 Sine, 3 Hydrostatic, 4 test
NumTerms=21;

if myCase == 1 %Uniform Load
    
    p= @(x,y) 1;
    p0=100000;
    
    La=.1; %Width
    Lb=.075; %Length
    
    %Material type and layup
    Mat_Types=[155E9,12.1E9,4.4E9,.248,-.018E-6,24.3E-6,.248,.458,24.3E-6];
        dlmwrite('materials.txt',Mat_Types);

        Angl=[30,0,0,30];
        Mtt=[1,1,1,1,];
        Thknss=[.001,.001,.001,.001,];


        for i=1:size(Angl,2)
            Big(i,1)=Mtt(i);
            Big(i,2)=Angl(i);
            Big(i,3)=Thknss(i);

        end
        dlmwrite('laminate.txt',Big);
        
        %Find the D matrix.
        Angles=Angl';
        Thicknesses=Thknss';
        Mat_Nums=Mtt';
        NL=size(Big,1);
        
        %use material numbers to create material property matrix
        for i=1:NL
            Mat_Props(i,:)=Mat_Types(Mat_Nums(i,1),:);
        end
        
        ABDx=ABD(Mat_Props,Angl',Thknss');
        D=ABDx(4:6,4:6);
        
        %Find the Fourier Series Coefficients.
        for n=1:2:NumTerms
            for m=1:2:NumTerms
                Qmn(m,n)=16*p0/(pi^2*m*n); %*(-1)^((m+n)/2-1);
                Wmn(m,n)=Qmn(m,n)/(D(1,1)*(m*pi/La)^4+2*(D(1,2)+2*D(3,3))*(m*pi/La)^2*(n*pi/Lb)^2+D(2,2)*(n*pi/Lb)^4);
            end
        end
elseif myCase == 2 %Sine Loading
    
    p= @(x,y) 100000;
    p0=100000;
    
    La=.1; %Width
    Lb=.075; %Length
    
    %Material type and layup
    Mat_Types=[155E9,12.1E9,4.4E9,.248,-.018E-6,24.3E-6,.248,.458,24.3E-6];
        dlmwrite('materials.txt',Mat_Types);

        Angl=[0,90,90,0];
        Mtt=[1,1,1,1,];
        Thknss=[.00015,.00015,.00015,.00015,];


        for i=1:size(Angl,2)
            Big(i,1)=Mtt(i);
            Big(i,2)=Angl(i);
            Big(i,3)=Thknss(i);

        end
        dlmwrite('laminate.txt',Big);
        
        %Find the D matrix.
        Angles=Angl';
        Thicknesses=Thknss';
        Mat_Nums=Mtt';
        NL=size(Big,1);
        
        %use material numbers to create material property matrix
        for i=1:NL
            Mat_Props(i,:)=Mat_Types(Mat_Nums(i,1),:);
        end
        
        ABDx=ABD(Mat_Props,Angl',Thknss');
        D=ABDx(4:6,4:6);
        
        %Find the Fourier Series Coefficients.
        Qmn(1,1)=p0;
        Qmn(2:NumTerms,2:NumTerms)=0;
        Wmn(1,1)=Qmn(1,1)/(D(1,1)*(1*pi/La)^4+2*(D(1,2)+2*D(3,3))*(1*pi/La)^2*(1*pi/Lb)^2+D(2,2)*(1*pi/Lb)^4);
        Wmn(2:NumTerms,2:NumTerms)=0;

elseif myCase == 3 % Hydrostatic loading        
        p= @(x,y) 100000;
    p0=100000;
    
    La=.15; %Width
    Lb=.05; %Length
    
    %Material type and layup
    Mat_Types=[155E9,12.1E9,4.4E9,.248,-.018E-6,24.3E-6,.248,.458,24.3E-6];
        dlmwrite('materials.txt',Mat_Types);

        Angl=[0,90,90,0];
        Mtt=[1,1,1,1,];
        Thknss=[.00015,.00015,.00015,.00015,];


        for i=1:size(Angl,2)
            Big(i,1)=Mtt(i);
            Big(i,2)=Angl(i);
            Big(i,3)=Thknss(i);

        end
        dlmwrite('laminate.txt',Big);
        
        %Find the D matrix.
        Angles=Angl';
        Thicknesses=Thknss';
        Mat_Nums=Mtt';
        NL=size(Big,1);
        
        %use material numbers to create material property matrix
        for i=1:NL
            Mat_Props(i,:)=Mat_Types(Mat_Nums(i,1),:);
        end
        
        ABDx=ABD(Mat_Props,Angl',Thknss');
        D=ABDx(4:6,4:6);
        
        %Find the Fourier Series Coefficients.
        for n=1:2:NumTerms
            for m=1:NumTerms
                Qmn(m,n)=8*p0/(pi^2*m*n)*(-1)^(m+1);
                Wmn(m,n)=Qmn(m,n)/(D(1,1)*(m*pi/La)^4+2*(D(1,2)+2*D(3,3))*(m*pi/La)^2*(n*pi/Lb)^2+D(2,2)*(n*pi/Lb)^4);
            end
        end
        
elseif myCase==4
    p= @(x,y) 1;
    p0=1;
    
    La=.5; %Width
    Lb=1; %Length
    a_want=10
    b_want=10
    
    %Material type and layup
    Mat_Types=[155E9,12.1E9,4.4E9,.248,-.018E-6,24.3E-6,.248,.458,24.3E-6];
        dlmwrite('materials.txt',Mat_Types);

        Angl=[90,0,0,90];
        Mtt=[1,1,1,1,];
        Thknss=[.00015,.00015,.00015,.00015,];


        for i=1:size(Angl,2)
            Big(i,1)=Mtt(i);
            Big(i,2)=Angl(i);
            Big(i,3)=Thknss(i);

        end
        dlmwrite('laminate.txt',Big);
        
        %Find the D matrix.
        Angles=Angl';
        Thicknesses=Thknss';
        Mat_Nums=Mtt';
        NL=size(Big,1);
        
        %use material numbers to create material property matrix
        for i=1:NL
            Mat_Props(i,:)=Mat_Types(Mat_Nums(i,1),:);
        end
        
        ABDx=ABD(Mat_Props,Angl',Thknss');
        D=ABDx(4:6,4:6);
        
        %Find the Fourier Series Coefficients.
        for n=1:2:NumTerms
            for m=1:2:NumTerms
                Qmn(m,n)=16*p0/(pi^2*m*n); %*(-1)^((m+n)/2-1);
                Wmn(m,n)=Qmn(m,n)/(D(1,1)*(m*pi/La)^4+2*(D(1,2)+2*D(3,3))*(m*pi/La)^2*(n*pi/Lb)^2+D(2,2)*(n*pi/Lb)^4);
            end
        end
end





%% Iterative Part
a_res=20; %number of nodes in the x
b_res=20;  %number of nodes in the y
for acount=1:a_res+1
    xPt=La*(acount-1)/a_res;
    for bcount=1:b_res+1
        yPt=Lb*(bcount-1)/b_res;
        
        
        %Find strains/curvatures based on location and load.
        Kap_x=0;
        Kap_y=0;
        Kap_xy=0;
        w=0;
        L=0;
        for mc=1:NumTerms
            for nc=1:NumTerms
                L=L+Qmn(mc,nc)*sin(mc*pi*xPt/La)*sin(nc*pi*yPt/Lb); %Calculate the Fourier approximation of the load
                w=w+(Wmn(mc,nc)*sin(mc*pi*xPt/La)*sin(nc*pi*yPt/Lb)); %  " " " of the force? 
                Kap_x=Kap_x+(Wmn(mc,nc)*(mc*pi/La)^2*sin(mc*pi*xPt/La)*sin(nc*pi*yPt/Lb)); % " " " of the Curvatures
                Kap_y=Kap_y+(Wmn(mc,nc)*(nc*pi/Lb)^2*sin(mc*pi*xPt/La)*sin(nc*pi*yPt/Lb));
                Kap_xy=Kap_xy+(-2)*(Wmn(mc,nc)*(mc*pi/La)*(nc*pi/Lb)*cos(mc*pi*xPt/La)*cos(nc*pi*yPt/Lb));
            end
        end
        
       %Find the most distant z location 
       zmax=-sum(Thknss)/2;
       
       %Since midplane strains are 0 Epsilon=Kappa*z
       Epsi_x=Kap_x*zmax;
       Epsi_y=Kap_y*zmax;
       Epsi_xy=Kap_xy*zmax;
       Epsi=[Epsi_x;Epsi_y;Epsi_xy];
       
        
        xsteps(acount)=xPt;
        ysteps(bcount)=yPt;
        eeK(1:3)=0;
        eeK(4:6)=[Kap_x,Kap_y,Kap_xy];
        dlmwrite('strains.txt', eeK) ;
        Proj4mod
        if (acount==a_want && bcount==b_want)
            bylayer(Mat_Props,Angles,Thicknesses,eeK)
            for r=1:8
                figure(r)
                axis auto
            end
        end
            
       Lxy(acount,bcount)=L;
       wxy(acount,bcount)=w;
       Kapx(acount,bcount)=Kap_x;
       Kapy(acount,bcount)=Kap_y;
       Kapxy(acount,bcount)=Kap_xy;
       
       
       sig_bot=off_axis_stress(:,3:5);
       sig_bot_on=on_axis_stress(:,3:5);
       
       sigbotx(acount,bcount,:)=sig_bot(:,1);
       sigboty(acount,bcount,:)=sig_bot(:,2);
       sigbotxy(acount,bcount,:)=sig_bot(:,3);
       
       sigbot1(acount,bcount,:)=sig_bot_on(:,1);
       sigbot2(acount,bcount,:)=sig_bot_on(:,2);
       sigbot12(acount,bcount,:)=sig_bot_on(:,3);
       
       eps1(acount,bcount,:)=on_axis_strain(:,3);
       eps2(acount,bcount,:)=on_axis_strain(:,4);
       eps12(acount,bcount,:)=on_axis_strain(:,5);
    end
end

%% Max Finding
[defmaxNum, defmaxIndex] = max(abs(wxy(:)));
[defx, defy] = ind2sub(size(wxy), defmaxIndex);

[kapxmaxNum, kapxmaxIndex] = max(abs(Kapx(:)));
[kapxx, kapxy] = ind2sub(size(Kapx), kapxmaxIndex);

[kapymaxNum, kapymaxIndex] = max(abs(Kapy(:)));
[kapyx, kapyy] = ind2sub(size(Kapy), kapymaxIndex);

[kapxymaxNum, kapxymaxIndex] = max(abs(Kapxy(:)));
[kapxyx, kapxyy] = ind2sub(size(Kapxy), kapxymaxIndex);

[sig1maxNum, sig1maxIndex] = max(abs(sigbot1(:)));
[sig1x, sig1y, sig1L] = ind2sub(size(sigbot1), sig1maxIndex);

[sig2maxNum, sig2maxIndex] = max(abs(sigbot2(:)));
[sig2x, sig2y, sig2L] = ind2sub(size(sigbot2), sig2maxIndex);

[sig12maxNum, sig12maxIndex] = max(abs(sigbot12(:)));
[sig12x, sig12y, sig12L] = ind2sub(size(sigbot12), sig12maxIndex);

[eps1maxNum, eps1maxIndex] = max(abs(eps1(:)));
[eps1x, eps1y, eps1L] = ind2sub(size(eps1), eps1maxIndex);

[eps2maxNum, eps2maxIndex] = max(abs(eps2(:)));
[eps2x, eps2y, eps2L] = ind2sub(size(eps2), eps2maxIndex);

[eps12maxNum, eps12maxIndex] = max(abs(eps12(:)));
[eps12x, eps12y, eps12L] = ind2sub(size(eps12), eps12maxIndex);

%% Plotting and outputs
%output the Qs
disp('Computer Assignment #5')
disp('Laminated Plates')
disp('Austin Stoker')
disp(' ')
disp('This program uses Fourier Series to calculate')
disp('laminate curvatures of a square plate')
disp('it passes these curvatures to another program')
disp('which returns the stress and strain')
disp('this program finds the maximum stresses, strains')
disp('and displacements and plots them for the layer')
disp('which contains the max value.')
disp(' ')
disp('========================================================')
disp(' ')
if myCase==1
    tmp='Uniform Loading';
elseif myCase==2
    tmp='Sinusoidal Loading';
elseif myCase==3
    tmp='Hydrostatic Loading';
end 
disp('This is for a square plate')
disp(tmp)
disp(strcat('p0 =',num2str(p0)))
disp(strcat('Width a = ',num2str(La)))
disp(strcat('Length b = ',num2str(Lb)))
disp(' ')
    
    
format short g
for g=1:NL
    Q_= Qbar(Mat_Props(g,:),Angl(g));
    temp=strcat('Qbar of Layer #',num2str(g));
    disp(temp)
    disp(Q_)
    disp(' ')
end
disp('========================================================')

disp(' ')
disp('D Matrix')
disp(D_mat)
disp(' ')


%           =========  The Load  ==============
di=max(La,Lb);

figure(50)
surf(xsteps,ysteps,Lxy')
title('The Fourier approximation of the Load')
axis([0 di 0 di])
view(3)


%           Displacements  =========================
%Plot
figure(1)
surf(xsteps,ysteps,wxy')
title('Displacements')
xlabel('x')
ylabel('y')
axis([0 di 0 di])
view(3)
%Text
Mat_Out=[xsteps(defx),ysteps(defy),wxy(defx,defy)];
disp('Max Deflection')
disp('         x            y         Deflection')
disp(Mat_Out)
disp(' ')


disp('=====================================================')
%=================== Kappa X  ==================

figure(2)
surf(xsteps,ysteps,Kapx')
title('Kappa x')
xlabel('x')
ylabel('y')
axis([0 di 0 di])
view(3)

Mat_Out=[xsteps(kapxx),ysteps(kapxy),kapxmaxNum];
disp('Kappa X')
disp('         x            y         Max Kappa X')
disp(Mat_Out)
disp(' ')

%  =================   Kappa y    =================
figure(3)
surf(xsteps,ysteps,Kapy')
title('Kappa Y')
xlabel('x')
ylabel('y')
axis([0 di 0 di])
view(3)

Mat_Out=[xsteps(kapyx),ysteps(kapyy),kapymaxNum];
disp('Kappa Y')
disp('         x            y         Max Kappa Y')
disp(Mat_Out)
disp(' ')


%=================== Kappa XY  ==================
figure(4)
surf(xsteps,ysteps,Kapxy')
title('Kappa xy')
xlabel('x')
ylabel('y')
axis([0 di 0 di])
view(3)

Mat_Out=[xsteps(kapyx),ysteps(kapyy),kapymaxNum];
disp('Kappa XY')
disp('         x            y         Max Kappa XY')
disp(Mat_Out)
disp(' ')

disp('=====================================================')

%=================== Epsilon 1 ===================
figure(5)
surf(xsteps,ysteps,eps1(:,:,eps1L)')
title('Epsilon 1')
xlabel('x')
ylabel('y')
axis([0 di 0 di])
view(3)
hold on
plot3(xsteps(eps1x),ysteps(eps1y),eps1(eps1x,eps1y,eps1L),'o','MarkerFaceColor','r','MarkerSize',15)
temp1= num2str(eps1(eps1x,eps1y,eps1L));
temp2= num2str(xsteps(eps1x));
temp3= num2str(ysteps(eps1y));
if cos(eps1L*pi)==-1
    temp=strcat('top of layer #',num2str((eps1L+1)/2));
else
    temp=strcat('bottom of layer #',num2str((eps1L)/2));
end
temp4=strcat(' (',temp2,',',temp3,',',temp1,')',temp);
text(xsteps(eps1x),ysteps(eps1y),eps1(eps1x,eps1y,eps1L),temp4,'EdgeColor','red','FontSize',12,'FontWeight','bold','BackgroundColor','w','HorizontalAlignment','right')


Mat_Out=[xsteps(eps1x),ysteps(eps1y),eps1maxNum];
disp('Max Epsilon 1')
disp(temp)
disp('         x            y         Epsilon 1')
disp(Mat_Out)
disp(' ')

%======================  Epsilon 2 ====================

figure(6)
surf(xsteps,ysteps,eps2(:,:,eps2L)')
title('Epsilon 2')
xlabel('x')
ylabel('y')
axis([0 di 0 di])
view(3)
hold on
plot3(xsteps(eps2x),ysteps(eps2y),eps2(eps2x,eps2y,eps2L),'o','MarkerFaceColor','r','MarkerSize',15)
temp1= num2str(eps2(eps2x,eps2y,eps2L));
temp2= num2str(xsteps(eps2x));
temp3= num2str(ysteps(eps2y));
if cos(eps2L*pi)==-1
    temp=strcat('top of layer #',num2str((eps2L+1)/2));
else
    temp=strcat('bottom of layer #',num2str((eps2L)/2));
end
temp4=strcat('(',temp2,',',temp3,',',temp1,') ',temp);
text(xsteps(eps2x),ysteps(eps2y),eps2(eps2x,eps2y,eps2L),temp4,'EdgeColor','red','FontSize',12,'FontWeight','bold','BackgroundColor','w','HorizontalAlignment','right')

Mat_Out=[xsteps(eps2x),ysteps(eps2y),eps2maxNum];
disp('Max Epsilon 2')
disp(temp)
disp('         x            y         Epsilon 2')
disp(Mat_Out)
disp(' ')


%====================== Gamma 12 ====================
figure(7)
surf(xsteps,ysteps,eps12(:,:,eps12L)')
title('Gamma 12')
xlabel('x')
ylabel('y')
axis([0 di 0 di])
view(3)
hold on
plot3(xsteps(eps12x),ysteps(eps12y),eps12(eps12x,eps12y,eps12L),'o','MarkerFaceColor','r','MarkerSize',15)
temp1= num2str(eps12(eps12x,eps12y,eps12L));
temp2= num2str(xsteps(eps12x));
temp3= num2str(ysteps(eps12y));
if cos(eps12L*pi)==-1
    temp=strcat('top of layer #',num2str((eps12L+1)/2));
else
    temp=strcat('bottom of layer #',num2str((eps12L)/2));
end
temp4=strcat('(',temp2,',',temp3,',',temp1,') ',temp);
text(xsteps(eps12x),ysteps(eps12y),eps12(eps12x,eps12y,eps12L),temp4,'EdgeColor','red','FontSize',12,'FontWeight','bold','BackgroundColor','w','HorizontalAlignment','right')
Mat_Out=[xsteps(eps12x),ysteps(eps12y),eps12maxNum];
disp('Max Gamma 12')
disp(temp)
disp('         x            y         Gamma 2')
disp(Mat_Out)
disp(' ')


disp('=====================================================')


%  ================ Sigma 1    =========================
%Plot
figure(8)
surf(xsteps,ysteps,sigbot1(:,:,sig1L)')
title('Sigma 1')
xlabel('x')
ylabel('y')
axis([0 di 0 di])
view(3)

%text
if cos(sig1L*pi)==-1
    temp=strcat('top of layer #',num2str((sig1L+1)/2));
else
    temp=strcat('bottom of layer #',num2str((sig1L)/2));
end
Mat_Out=[xsteps(sig1x),ysteps(sig1y),sig1maxNum];
disp('Max Sigma 1')
disp(temp)
disp('         x            y         Sigma 1')
disp(Mat_Out)
disp(' ')


%=============    Sigma 2 ====================
%plot
figure(9)
surf(xsteps,ysteps,sigbot2(:,:,sig2L)')
title('Sigma 2')
xlabel('x')
ylabel('y')
axis([0 di 0 di])
view(3)

%text
if cos(sig2L*pi)==-1
    temp=strcat('top of layer #',num2str((sig2L+1)/2));
else
    temp=strcat('bottom of layer #',num2str((sig2L)/2));
end
Mat_Out=[xsteps(sig2x),ysteps(sig2y),sig2maxNum];
disp('Max Sigma 2')
disp(temp)
disp('         x            y         Sigma 2')
disp(Mat_Out)
disp(' ')



%==================== Tau 12 ====================
%Plot
figure(10)
surf(xsteps,ysteps,sigbot12(:,:,1)')
title('Tau 12')
xlabel('x')
ylabel('y')
axis([0 di 0 di])
view(3)

%text
if cos(sig12L*pi)==-1
    temp=strcat('top of layer #',num2str((sig12L+1)/2));
else
    temp=strcat('bottom of layer #',num2str((sig12L)/2));
end
Mat_Out=[xsteps(sig12x),ysteps(sig12y),sig12maxNum];
disp('Max Tau 12')
disp(temp)
disp('         x            y         Tau 12')
disp(Mat_Out)
disp(' ')



disp('===================================================')


