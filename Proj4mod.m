%%Austin Stoker

%Composites Project #4

%Take material properties and lamina orientations from file
%as well as strains and curvatures or loads and moments

%Print Q,S,Qbar,Sbar,ABD,abd matrices
%Smeared properties
%and Strains and curvatures or loads and moments whichever wasn't given.

%% Setup



%============================================
%============================================
close all
clear all
clc

%get the Material types list
Mat_Types=dlmread('materials.txt');

%get the layout of the laminate
Big=dlmread('laminate.txt');

%Find how many layers
NL=size(Big,1);

%Find how many material types
NT=size(Mat_Types,1);


%display the material types list
for i=1:size(Mat_Types,1)
    MatList(i,1)=i;
    MatList(i,2:5)=Mat_Types(i,1:4);
end


%get angles thicknesses and material numbers
Angles=Big(:,2);
Thicknesses=Big(:,3);
Mat_Nums=Big(:,1);
H=sum(Thicknesses);

%use material numbers to create material property matrix
for i=1:NL
    Mat_Props(i,:)=Mat_Types(Mat_Nums(i,1),:);
end

dT=-200;
dM=0;



%set up the layout of the composite
for i=1:size(Big,1)
    layout(i,1)=i;
    layout(i,2:4)=Big(i,1:3);
end



%Check for Loads.txt
%if it exists use the loads to find strains and curvatures
LStatus=0;
SStatus=0;
Strains=0;
Stresses=0;
LStatus=exist('loads.txt','file');
    if LStatus==2 
        %loads are given so display them
        Loads=dlmread('loads.txt');
        N_x=Loads(1);
        N_y=Loads(2);
        N_xy=Loads(3);
        M_x=Loads(4);
        M_y=Loads(5);
        M_xy=Loads(6);
        
    else
        SStatus=exist('strains.txt','file');
        if SStatus==2 
            %Strains are given so load and display them
            Strains=dlmread('strains.txt');
            eps_x=Strains(1);
            eps_y=Strains(2);
            eps_xy=Strains(3);
            Kappa_x=Strains(4);
            Kappa_y=Strains(5);
            Kappa_xy=Strains(6);
            
            
        end
    end
    

    



%% Thermal and Hygral effects
 % Read in the thermal and hygral properties
alpha_12=Mat_Props(:,5:6);
alpha_12(:,3)=0;
beta_12=Mat_Props(:,7:8);
beta_12(:,3)=0;


Tt=0;
for k=1:NL
    Tt=Thicknesses(k) + Tt;
end

z(1)=0-Tt/2;
for k=1:NL
   z(k+1)=Thicknesses(k) + z(k);
end

N_hat_T=0;
for k=1:NL
    alpha_xy(k,:)=(T(Angles(k))*alpha_12(k,1:3)')';
    alpha_xy(k,3)=-2*alpha_xy(k,3);
    N_hat_T=Qbar(Mat_Props(k,1:4),Angles(k))*alpha_xy(k,:)'*Thicknesses(k)+ N_hat_T;
end

N_hat_M=0;
for k=1:NL
    beta_xy(k,:)=(T(Angles(k))*beta_12(k,1:3)')';
    beta_xy(k,3)=-2*beta_xy(k,3);
    N_hat_M=Qbar(Mat_Props(k,1:4),Angles(k))*beta_xy(k,:)'*Thicknesses(k)+ N_hat_M;
end

M_hat_T=0;
for k=1:NL
    alpha_xy(k,:)=(T(Angles(k))*alpha_12(k,1:3)')';
    alpha_xy(k,3)=-2*alpha_xy(k,3);
    M_hat_T=Qbar(Mat_Props(k,1:4),Angles(k))*alpha_xy(k,:)'*((z(k+1))^2 - (z(k))^2)/2+ M_hat_T;
end

M_hat_M=0;
for k=1:NL
    beta_xy(k,:)=(T(Angles(k))*beta_12(k,1:3)')';
    beta_xy(k,3)=-2*beta_xy(k,3);
    M_hat_M=Qbar(Mat_Props(k,1:4),Angles(k))*beta_xy(k,:)'*((z(k+1))^2 - (z(k))^2)/2+ M_hat_M;
end






%% The ABD
%Get and display the ABD matrix
ABD_Mat=ABD(Mat_Props,Angles,Thicknesses);
A_mat=ABD_Mat(1:3,1:3);
B_mat=ABD_Mat(1:3,4:6);
D_mat=ABD_Mat(4:6,4:6);




%Get and display the abd matrix
abd_Mat=ABD_Mat^-1;
a_mat=abd_Mat(1:3,1:3);
b_mat=abd_Mat(1:3,4:6);
d_mat=abd_Mat(4:6,4:6);



%% Find Loads or Moments
%if the Loads were give find the midplane Strains
if LStatus==2
    AppliedLoads=Loads';      

    ThermLoads(1)=N_hat_T(1)*dT+N_hat_M(1)*dM;
    ThermLoads(2)=N_hat_T(2)*dT+N_hat_M(2)*dM;
    ThermLoads(3)=N_hat_T(3)*dT+N_hat_M(3)*dM;
    ThermLoads(4)=M_hat_T(1)*dT+M_hat_M(1)*dM;
    ThermLoads(5)=M_hat_T(2)*dT+M_hat_M(2)*dM;
    ThermLoads(6)=M_hat_T(3)*dT+M_hat_M(3)*dM;

    TotalLoads=AppliedLoads+ThermLoads';
    

    Strains=abd_Mat*TotalLoads;
    eps_x=Strains(1);
    eps_y=Strains(2);
    eps_xy=Strains(3);
    Kappa_x=Strains(4);
    Kappa_y=Strains(5);
    Kappa_xy=Strains(6);

   
end
% 
% 
% %if the Strains were given, calculate the Loads
if SStatus==2
    Loads=ABD_Mat*Strains';
    N_x=Loads(1);
    N_y=Loads(2);
    N_xy=Loads(3);
    M_x=Loads(4);
    M_y=Loads(5);
    M_xy=Loads(6);

    
AppliedLoads=Loads;      

ThermLoads(1)=N_hat_T(1)*dT+N_hat_M(1)*dM;
ThermLoads(2)=N_hat_T(2)*dT+N_hat_M(2)*dM;
ThermLoads(3)=N_hat_T(3)*dT+N_hat_M(3)*dM;
ThermLoads(4)=M_hat_T(1)*dT+M_hat_M(1)*dM;
ThermLoads(5)=M_hat_T(2)*dT+M_hat_M(2)*dM;
ThermLoads(6)=M_hat_T(3)*dT+M_hat_M(3)*dM;

TotalLoads=AppliedLoads+ThermLoads';
end




%% Smeared Properties

%Calculate the smeared properties
  

E_x=1/(a_mat(1,1)*H);
E_y=1/(a_mat(2,2)*H);
G_xy=1/(a_mat(3,3)*H);

v_xy=-a_mat(1,2)/a_mat(1,1);
v_yx=-a_mat(1,2)/a_mat(2,2);


n_x_xy=a_mat(1,3)/a_mat(3,3);
n_y_xy=a_mat(2,3)/a_mat(3,3);

n_xy_x=a_mat(1,3)/a_mat(1,1);
n_xy_y=a_mat(2,3)/a_mat(2,2);

alpha_bar=a_mat*N_hat_T;
alpha_bar_x=alpha_bar(1);
alpha_bar_y=alpha_bar(2);
alpha_bar_xy=alpha_bar(3);

%% Project #3 begins
% Calculating the various stresses based on our midplane strains and
% curvatures

eK=Strains;
NL=size(Angles,1); %number of layers



%Total thickness of the laminate
Tt=0;
for k=1:NL
    Tt=Thicknesses(k) + Tt;
end


%positions of lamina edges, z(0) in the book = z(1) here because
%matlab doesn't do 0 index
z(1)=0-Tt/2;
for k=1:NL
   z(k+1)=Thicknesses(k) + z(k);
end

%make a list of z's that includes 2 points for all interior z's
z2(1)=z(1);
Angles2(1)=Angles(1);
Angles2(2)=Angles(1);
Mat_Props2(1,:)=Mat_Props(1,:);
Mat_Props2(2,:)=Mat_Props(1,:);
alpha_xy2(1,:)=alpha_xy(1,:);
alpha_xy2(2,:)=alpha_xy(1,:);

for i=2:NL
    z2(2*i-2)=z(i);
    z2(2*i-1)=z(i);
    Angles2(2*i-1)=Angles(i);
    Angles2(2*i)=Angles(i);
    Mat_Props2(2*i-1,:)=Mat_Props(i,:);
    Mat_Props2(2*i,:)=Mat_Props(i,:);
    alpha_xy2(2*i-1,:)=alpha_xy(i,:);
    alpha_xy2(2*i,:)=alpha_xy(i,:);
end
z2(2*NL)=z(NL+1);



%Find the xy Strains
for k=1:2*NL
Mech_Eps_x(k) = eK(1)+ z2(k)*eK(4)-alpha_xy2(k,1)*dT;
Mech_Eps_y(k) = eK(2)+ z2(k)*eK(5)-alpha_xy2(k,2)*dT;
Mech_Tau_xy(k) = eK(3)+ z2(k)*eK(6)-alpha_xy2(k,3)*dT;

Meas_Eps_x(k) = eK(1)+ z2(k)*eK(4);
Meas_Eps_y(k) = eK(2)+ z2(k)*eK(5);
Meas_Tau_xy(k) = eK(3)+ z2(k)*eK(6);
end


%Convert to 1-2 Strains
Mech_Strains_xy=[Mech_Eps_x(1:end);Mech_Eps_y(1:end);Mech_Tau_xy(1:end)];
Meas_Strains_xy=[Meas_Eps_x(1:end);Meas_Eps_y(1:end);Meas_Tau_xy(1:end)];

R = [1,0,0;0,1,0;0,0,2];
for i=1:2*NL
    Mech_Strains_12(:,i)= R*T(Angles2(i))*R'*Mech_Strains_xy(:,i);
    Meas_Strains_12(:,i)= R*T(Angles2(i))*R'*Meas_Strains_xy(:,i);
end

%Find xy Stresses
for i=1:2*NL
    Stresses_xy(:,i)=Qbar(Mat_Props2(i,1:4),Angles2(i))*Mech_Strains_xy(:,i);
end


%Convert to 1-2 Stresses
for i=1:2*NL
    Stresses_12(:,i)= T(Angles2(i))*Stresses_xy(:,i);
end

Sig_x=Stresses_xy(1,:);
Sig_y=Stresses_xy(2,:);
Tau_xy=Stresses_xy(3,:);
Sig_1=Stresses_12(1,:);
Sig_2=Stresses_12(2,:);
Tau_12=Stresses_12(3,:);

Meas_Eps_x=Meas_Strains_xy(1,:);
Meas_Eps_y=Meas_Strains_xy(2,:);
Meas_Gam_xy=Meas_Strains_xy(3,:);
Meas_Eps_1=Meas_Strains_12(1,:);
Meas_Eps_2=Meas_Strains_12(2,:);
Meas_Gam_12=Meas_Strains_12(3,:);

%=========================================
%Get rid of little tiny values

max_stress=max(max(max(abs(Stresses_xy),abs(Stresses_12))));
for i=1:3;
   for j=1:2*NL;
       if abs(Stresses_xy(i,j))< abs(max_stress)*10E-10
           Stresses_xy(i,j)=0;
       end
       if abs(Stresses_12(i,j))< abs(max_stress)*10E-10
           Stresses_12(i,j)=0;
       end
   end
end

% max_strain=max(max(max(abs(Meas_Strains_xy),abs(Meas_Strains_12))));
% for i=1:3;
%    for j=1:2*NL;
%        if abs(Meas_Strains_xy(i,j))< abs(max_strain)*10E-10
%            Meas_Strains_xy(i,j)=0;
%        end
%        if abs(Meas_Strains_12(i,j))< abs(max_strain)*10E-10
%            Meas_Strains_12(i,j)=0;
%        end
%    end
% end


%Smash it all together into a matrix
for i=1:NL
    %put in the layer # for the top
    off_axis_stress(2*i-1,1)=i;
    %put a 1 in so 'Top' will be displayed
    off_axis_stress(2*i-1,2)=1;
    %put in the layer # for the bottom
    off_axis_stress(2*i,1)=i;
    %put a 1 in so 'Bottom' will be displayed
    off_axis_stress(2*i,2)=-1;
    
    off_axis_strain(2*i-1,1)=i;
    off_axis_strain(2*i-1,2)=1;
    off_axis_strain(2*i,1)=i;
    off_axis_strain(2*i,2)=-1;
    
    on_axis_stress(2*i-1,1)=i;
    on_axis_stress(2*i-1,2)=1;
    on_axis_stress(2*i,1)=i;
    on_axis_stress(2*i,2)=-1;
    
    on_axis_strain(2*i-1,1)=i;
    on_axis_strain(2*i-1,2)=1;
    on_axis_strain(2*i,1)=i;
    on_axis_strain(2*i,2)=-1;
end
    off_axis_stress(:,3:5)= Stresses_xy';
    off_axis_strain(:,3:5)= Meas_Strains_xy';
    
    on_axis_stress(:,3:5)= Stresses_12';
    on_axis_strain(:,3:5)= Meas_Strains_12';
   
    
%% Failure Tests    
    
    fail_crit=dlmread('fail_crit.txt');
    
    failure(:,1:2)=on_axis_stress(:,1:2);
   
    F1=(1/fail_crit(1)+1/(fail_crit(2)));
    F11=-1/(fail_crit(1)*fail_crit(2));
    
    F2=(1/fail_crit(3)+1/(fail_crit(4)));
    F22=-1/(fail_crit(3)*fail_crit(4));
        
    F66=1/(fail_crit(5)^2);
        
    for i=1:2*NL
        for j=3:4
        if on_axis_stress(i,j)<0 % compression
            failure(i,j)=-on_axis_stress(i,j)/fail_crit(2*(j-2));
        else %Tension or no load
            failure(i,j)=on_axis_stress(i,j)/fail_crit(2*(j-2)-1);
        end
        end
        failure(i,5)=abs(on_axis_stress(i,5))/fail_crit(5);
        %Tsai Wu
        
        
        failure(i,6)=F1*on_axis_stress(i,3)+F2*on_axis_stress(i,4)+F11*(on_axis_stress(i,3)^2)+F22*(on_axis_stress(i,4)^2)+F66*(on_axis_stress(i,5)^2)-sqrt((F11*F22))*on_axis_stress(i,3)*on_axis_stress(i,4);
    end    
%% Through The Thickness    
% Find the eps_bar_z
eps_bar_z=0;
for i=1:NL
    
    eps_bar_z=(1/H)*((-Mat_Props(i,7)/Mat_Props(i,1))*Stresses_12(1,2*i-1)- (Mat_Props(i,8)/Mat_Props(i,2))*Stresses_12(2,2*i-1))*Thicknesses(i) + eps_bar_z;
    

end

eps_bar_z_T=0;
for i=1:NL
    
    eps_bar_z_T=(1/H)*(Mat_Props(i,9)*dT+(-Mat_Props(i,7)/Mat_Props(i,1))*Stresses_12(1,2*i-1)- (Mat_Props(i,8)/Mat_Props(i,2))*Stresses_12(2,2*i-1))*Thicknesses(i) + eps_bar_z_T;
    

end


v_xz=0;
v_yz=0;
alpha_bar_z=0;


for i=1:NL
    S13 = -Mat_Props(i,7)/Mat_Props(i,1);
    S23 = -Mat_Props(i,8)/Mat_Props(i,2);
    Q_Mat = Qreduced(Mat_Props(i,1),Mat_Props(i,2),Mat_Props(i,3),Mat_Props(i,4));
    

    Q11=Q_Mat(1,1);
    Q12 = Q_Mat(1,2);
    Q22 = Q_Mat(2,2);
    a11=a_mat(1,1);
    a12=a_mat(1,2);
    a22=a_mat(2,2);
    a16=a_mat(1,3);
    m=cos(Angles(i)*pi/180);
    n=sin(Angles(i)*pi/180);
    alpha_3=Mat_Props(i,9);
    
  
    
    v_xz=(-1/(a_mat(1,1)*H))*((S13*(Q11*m^2+Q12*n^2)+S23*(Q12*m^2+Q22*n^2))*a11+ (S13*(Q11*n^2+Q12*m^2)+S23*(Q12*n^2+Q22*m^2))*a12+ (S13*(Q11-Q12)+S23*(Q12-Q22))*m*n*a16)*Thicknesses(i)  + v_xz;
    v_yz=(-1/(a_mat(2,2)*H))*((S13*(Q11*m^2+Q12*n^2)+S23*(Q12*m^2+Q22*n^2))*a12+ (S13*(Q11*n^2+Q12*m^2)+S23*(Q12*n^2+Q22*m^2))*a22+ (S13*(Q11-Q12)+S23*(Q12-Q22))*m*n*a16)*Thicknesses(i)  + v_yz;
    alpha_bar_z = alpha_bar_z + ((S13*(Q11*m^2+Q12*n^2)+S23*(Q12*m^2+Q22*n^2))*alpha_bar_x+ (S13*(Q11*n^2+Q12*m^2)+S23*(Q12*n^2+Q22*m^2))*alpha_bar_y + alpha_3 - (S13*(Q11*alpha_12(i,1)+Q12*alpha_12(i,2))+S23*(Q12*alpha_12(i,1)+Q22*alpha_12(i,2))));
end
    alpha_bar_z=  alpha_bar_z/NL;

    
