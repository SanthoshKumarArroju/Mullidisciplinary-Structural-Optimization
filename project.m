clear all
close all
clc
   %% FEA Analysis %%

       E = 73e9, r=2780, g=9.81 
       n1=[0,0];
       n2=[0.8,0];
       n3=[1.6,0];
       n4=[2.4,0];
       n5=[3.2,0];
       n6=[4,0];
       n7=[4,0.8];
       n8=[2.4,0.64];
       n9=[0.8,0.48];
       n10=[0,0.4];
       
       L1 = sqrt(((n2(1)-n1(1))^2+(n2(2)-n1(2)^2)));
       L2 = sqrt(((n3(1)-n2(1))^2+(n3(2)-n2(2)^2)));
       L3 = sqrt(((n4(1)-n3(1))^2+(n4(2)-n3(2)^2)));
       L4 = sqrt(((n5(1)-n4(1))^2+(n5(2)-n4(2)^2)));
       L5 = sqrt(((n6(1)-n5(1))^2+(n6(2)-n5(2)^2)));
       L6 = sqrt(((n7(1)-n6(1))^2+(n7(2)-n6(2)^2)));
       L7 = sqrt(((n8(1)-n7(1))^2+(n8(2)-n7(2)^2)));
       L8 = sqrt(((n9(1)-n8(1))^2+(n9(2)-n8(2)^2)));
       L9 = sqrt(((n10(1)-n9(1))^2+(n10(2)-n9(2)^2)));
       L10 = sqrt(((n10(1)-n1(1))^2+(n10(2)-n1(2)^2)));
       L11 = sqrt(((n9(1)-n1(1))^2+(n9(2)-n1(2)^2)));
       L12 = sqrt(((n9(1)-n2(1))^2+(n9(2)-n2(2)^2)));
       L13 = sqrt(((n9(1)-n3(1))^2+(n9(2)-n3(2)^2)));
       L14 = sqrt(((n3(1)-n8(1))^2+(n3(2)-n8(2)^2)));
       L15 = sqrt(((n8(1)-n4(1))^2+(n8(2)-n4(2)^2)));
       L16 = sqrt(((n8(1)-n5(1))^2+(n8(2)-n5(2)^2)));
       L17 = sqrt(((n7(1)-n5(1))^2+(n7(2)-n5(2)^2)));
       tetha1 = 0;
       tetha2 = 0;
       tetha3 = 0;
       tetha4 = 0;
       tetha5 = 0;
       tetha6 = 90;
       tetha7 = 180 + atand(0.16/(L4+L5));
       tetha8 = 180 + atand(0.16/(L2+L3));
       tetha9 = 180 + atand(0.08/L1);
       tetha10 = 0;
       tetha11 = atand(L12/L1);
       tetha12 = 0;
       tetha13 = 180 - atand(L12/L2);
       tetha14 = atand(L15/L3);
       tetha15 = 0;
       tetha16 = 180 - atand(L15/L4);
       tetha17 = atand(L6/L5);
      
       L = [L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12 L13 L14 L15 L16 L17]
       thet = [tetha1 tetha2 tetha3 tetha4 tetha5 tetha6 tetha7 tetha8 tetha9 tetha10 tetha11 tetha12 tetha13 tetha14 tetha15 tetha16 tetha17];
       R = [0.020 0.020 0.020 0.020 0.020 0.015 0.024 0.024 0.024 0.015 0.015 0.015 0.015 0.015 0.015 0.015 0.015]
%        A = (22/7)*R.^2
A = [-15788075667.9907,-11506646463.7620,-8078996746.54495,-5410180780.93450,-3405252831.52555,-1969267162.91298,-1007278039.69170,-424339726.456621,-125506487.802635,-15832588.3246507,-372292.617570435,15820134.7237020,127690429.104263,430184325.929210,1018247560.60364,1986825868.53265,3430864985.12133];
       I = (22/7)*R.^4 /4
       
       for i=1:17;  % Beam elements %
            k_localbeam(:,:,i) = [(E*A(i))/L(i),                   0,                 0, -(E*A(i))/L(i),                   0,                   0;
                                                                 0,  (12*E*I(i))/L(i)^3, (6*E*I(i))/L(i)^2,                                      0, -(12*E*I(i))/L(i)^3,   (6*E*I(i))/L(i)^2;
                                                                 0,   (6*E*I(i))/L(i)^2,   (4*E*I(i))/L(i),                                      0,  -(6*E*I(i))/L(i)^2,     (2*E*I(i))/L(i);
                            -(E*A(i))/L(i),                   0,                  0, (E*A(i))/L(i),                   0,                   0;
                                                                 0, -(12*E*I(i))/L(i)^3, -(6*E*I(i))/L(i)^2,                                      0,  (12*E*I(i))/L(i)^3, -(6*E*I(i))/L(i)^2;
                                                                 0,   (6*E*I(i))/L(i)^2,    (2*E*I(i))/L(i),                                      0,  -(6*E*I(i))/L(i)^2,    (4*E*I(i))/L(i)];


             M_localbeam(:,:,i)=(r*A(i)*L(i))*[2/6 0 0 1/6 0 0;
                                 0 156/420 (22*L(i))/420 0 54/420 (-13*L(i))/420;
                            0 (22*L(i))/420 (4*L(i)*L(i))/420 0 (13*L(i))/420 (-3*L(i)*L(i))/420;
                                                                1/6 0 0 2/6 0 0;
                                       0 54/420 (13*L(i))/420 0 156/420 (-22*L(i))/420;
                            0 (-13*L(i))/420 (-3*L(i)*L(i))/420 0 (-22*L(i))/420 (4*L(i)*L(i))/420];
                                                             
            dk_da_localbeam(:,:,i) = [(E)/L(i),                   0,                 0, -(E)/L(i),                   0,                   0;
                                                                 0,  (12*E*I(i))/L(i)^3, (6*E*I(i))/L(i)^2,                                      0, -(12*E*I(i))/L(i)^3,   (6*E*I(i))/L(i)^2;
                                                                 0,   (6*E*I(i))/L(i)^2,   (4*E*I(i))/L(i),                                      0,  -(6*E*I(i))/L(i)^2,     (2*E*I(i))/L(i);
                                    -(E)/L(i),                   0,                  0, (E)/L(i),                   0,                   0;
                                                                 0, -(12*E*I(i))/L(i)^3, -(6*E*I(i))/L(i)^2,                                      0,  (12*E*I(i))/L(i)^3, -(6*E*I(i))/L(i)^2;
                                                                 0,   (6*E*I(i))/L(i)^2,    (2*E*I(i))/L(i),                                      0,  -(6*E*I(i))/L(i)^2,    (4*E*I(i))/L(i)];

      dM_da_localbeam(:,:,i)= (r*L(i))*[2/6 0 0 1/6 0 0;
                                 0 156/420 (22*L(i))/420 0 54/420 (-13*L(i))/420;
                            0 (22*L(i))/420 (4*L(i)*L(i))/420 0 (13*L(i))/420 (-3*L(i)*L(i))/420;
                                                                1/6 0 0 2/6 0 0;
                                       0 54/420 (13*L(i))/420 0 156/420 (-22*L(i))/420;
                            0 (-13*L(i))/420 (-3*L(i)*L(i))/420 0 (-22*L(i))/420 (4*L(i)*L(i))/420];
                                                             
       T_beam(:,:,i)=[cosd(thet(i)) sind(thet(i)) 0 0 0 0;
                     -sind(thet(i)) cosd(thet(i)) 0 0 0 0;
                                              0 0 1 0 0 0;
                      0 0 0 cosd(thet(i)) sind(thet(i)) 0;
                     0 0 0 -sind(thet(i)) cosd(thet(i)) 0;
                                              0 0 0 0 0 1];
       k_beam(:,:,i) = T_beam(:,:,i)'*k_localbeam(:,:,i)*T_beam(:,:,i);    % k global for Beam elements %
       M_beam(:,:,i) = T_beam(:,:,i)'*M_localbeam(:,:,i)*T_beam(:,:,i);
       dk_da_beam(:,:,i) = T_beam(:,:,i)'*dk_da_localbeam(:,:,i)*T_beam(:,:,i) ;   % dk/dt  global Beam elements %    
       dM_da_beam(:,:,i) = T_beam(:,:,i)'*dM_da_localbeam(:,:,i)*T_beam(:,:,i) ;     
 
       end
       
   %% Global Stiffness Matrix (kg_global = kg_truss + kg_beam) %% 
      format short
      kg_global=zeros(30,30,1); dkg_global = zeros(30,30,17);
      Mg_global=zeros(30,30,1); dMg_global = zeros(30,30,17);
for i = 1:17;
     Dof = zeros(17,17); 
Dof = [1 2 3 4 5 6;4 5 6 7 8 9; 7 8 9 10 11 12; 10 11 12 13 14 15; 13 14 15 16 17 18; 16 17 18 19 20 21; 19 20 21 22 23 24; 22 23 24 25 26 27; 25 26 27 28 29 30;28 29 30 1 2 3;1 2 3 25 26 27; 4 5 6 25 26 27; 25 26 27 7 8 9; 7 8 9 22 23 24; 10 11 12 22 23 24; 22 23 24 13 14 15; 13 14 15 19 20 21];     
     kg_global(Dof(i,:),Dof(i,:),1)= kg_global(Dof(i,:),Dof(i,:),1)+ k_beam(:,:,i);
     dkg_global(Dof(i,:),Dof(i,:),i)= dkg_global(Dof(i,:),Dof(i,:),i)+ dk_da_beam(:,:,i);
     Mg_global(Dof(i,:),Dof(i,:),1)= Mg_global(Dof(i,:),Dof(i,:),1)+M_beam(:,:,i);
     dMg_global(Dof(i,:),Dof(i,:),i)= dMg_global(Dof(i,:),Dof(i,:),i)+dM_da_beam(:,:,i);

end 
     kg_beforeBC = kg_global;
     Mg_beforeBC = Mg_global;
     kg_beforeBC, Mg_beforeBC     
   %% Force Vector (Fg) %%
     Fg=zeros(30,1); Fg(14)= -7000; Fg(23) = -1360;     % considering Engine concentrated loads at 14, fuel weight at 23 %
  
   %% Apply Boundary conditions %% U[11,12,13,14]=0 
    
BC = 16; kg_global(BC,:)=0; kg_global(:,BC)=0; kg_global(BC,BC)=1; Fg(BC)=0;
BC = 17; kg_global(BC,:)=0; kg_global(:,BC)=0; kg_global(BC,BC)=1; Fg(BC)=0;
BC = 18; kg_global(BC,:)=0; kg_global(:,BC)=0; kg_global(BC,BC)=1; Fg(BC)=0;
BC = 19; kg_global(BC,:)=0; kg_global(:,BC)=0; kg_global(BC,BC)=1; Fg(BC)=0;
BC = 20; kg_global(BC,:)=0; kg_global(:,BC)=0; kg_global(BC,BC)=1; Fg(BC)=0;

BC = 16; Mg_global(BC,:)=0; Mg_global(:,BC)=0; Mg_global(BC,BC)=1; Fg(BC)=0;
BC = 17; Mg_global(BC,:)=0; Mg_global(:,BC)=0; Mg_global(BC,BC)=1; Fg(BC)=0;
BC = 18; Mg_global(BC,:)=0; Mg_global(:,BC)=0; Mg_global(BC,BC)=1; Fg(BC)=0;
BC = 19; Mg_global(BC,:)=0; Mg_global(:,BC)=0; Mg_global(BC,BC)=1; Fg(BC)=0;
BC = 20; Mg_global(BC,:)=0; Mg_global(:,BC)=0; Mg_global(BC,BC)=1; Fg(BC)=0;


     kg_afterBC = kg_global;
     Mg_afterBC = Mg_global;
     kg_afterBC, Mg_afterBC

 %% Displacement
    format short
    Ug = inv(kg_afterBC)*Fg
    
 %% Eigen Values and Eigen Vectors
  [V,D]=eig(kg_afterBC,Mg_afterBC)        % V: eigen vectors, D: eigen values %
  
 %% Natural Frequency and Mode shapes
  Natural_frequency= sqrt(D)/(2*pi)
  Mode_shapes = V  
 %% at the initila design point the constraint and weight function value A(i) are
 frequency(1) = 0; % initial freq is zero, free vibration
 weight = r*g*(A(1)*L(1)+A(2)*L(2)+A(3)*L(3)+A(4)*L(4)+A(5)*L(5)+A(6)*L(6)+A(7)*L(7)+A(8)*L(8)+A(9)*L(9)+A(10)*L(10)+A(11)*L(11)+A(12)*L(12)+A(13)*L(13)+A(14)*L(14)+A(15)*L(15)+A(16)*L(16)+A(17)*L(17))
 
 g1 = 1- Ug(2)/20e-03
 g2 = 1- Ug(5)/20e-03
 g3 = 1- Ug(8)/20e-03
 g4 = 1- Ug(11)/20e-03
 g5 = 1- Ug(14)/20e-03
 g6 = 1- Ug(23)/20e-03
 g7 = 1- frequency(1)/20e-03
  
  %% Sensitivity Analysis of weight of structure w.r.t. Design Variable 'R: radius'

 W = r*g*A.*L(i);
 for i=1:17
 Weight_Sensitivity(i) =r*g*(22/7)*2*R(i)*L(i);
 end
 Weight_Sensitivity
 
%% Sensitivity Analysis

%Constraint_1; g1 = Ug(2) < -20e-03 m (tip displacement Dof_2)
for i=1:17
    dUg_dU2 = [ 0;1 ;zeros(28,1)];  % partial diff g1 w.r.t Ug,  dUg/dUg(2)    
    dg1_da(:,:,i) = -(-1/20e-03)*(dUg_dU2)'*(inv(kg_afterBC))*(-dkg_global(:,:,i)*Ug);
end 
dg1_da
%Constraint_2; g2 = Ug(5) < -20e-03 m (displacement Dof_5)
for i=1:17
    dUg_dU5 = [ 0;0;0;0;0;1 ;zeros(24,1)];  % partial diff g1 w.r.t Ug,  dUg/dUg(6)    
    dg2_da(:,:,i) = -(-1/20e-03)*(dUg_dU5)'*(inv(kg_afterBC))*(-dkg_global(:,:,i)*Ug);
end 
dg2_da
%Constraint_3; g3 = Ug(8) < -20e-03 m (displacement Dof_8)
for i=1:17
    dUg_dU8 = [ zeros(7,1);1 ;zeros(22,1)];  % partial diff g1 w.r.t Ug,  dUg/dUg(11)    
    dg3_da(:,:,i) = -(-1/20e-03)*(dUg_dU8)'*(inv(kg_afterBC))*(-dkg_global(:,:,i)*Ug);
end 
dg3_da
%Constraint_4; g4 = Ug(11) < -20e-03 m (displacement Dof_11)
for i=1:17
    dUg_dU11 = [ zeros(10,1);1 ;zeros(19,1)];  % partial diff g1 w.r.t Ug,  dUg/dUg(1)    
    dg4_da(:,:,i) = -(-1/20e-03)*(dUg_dU11)'*(inv(kg_afterBC))*(-dkg_global(:,:,i)*Ug);
end 
dg4_da
%Constraint_5; g5 = Ug(14) < -10e-03 m (displacement Dof_14) 
for i=1:17
    dUg_dU14 = [ zeros(13,1);1 ;zeros(16,1)];  % partial diff g1 w.r.t Ug,  dUg/dUg(1)    
    dg5_da(:,:,i) = -(-1/10e-03)*(dUg_dU14)'*(inv(kg_afterBC))*(-dkg_global(:,:,i)*Ug);
end 
dg5_da
%Constraint_6; g6 = Ug(23) < -10e-03 m (displacement Dof_23)
for i=1:17
    dUg_dU23 = [zeros(22,1);1 ;zeros(7,1)];  % partial diff g1 w.r.t Ug,  dUg/dUg(23)    
    dg6_da(:,:,i) = -(-1/10e-03)*(dUg_dU23)'*(inv(kg_afterBC))*(-dkg_global(:,:,i)*Ug);
end 
dg6_da
%Constraint_7; g7 = freq(1) > 20 Hz (1st Natural Frequency > 20 Hz)
% V = eigen_vectors
df_dD1 = 1/(4*pi*sqrt(D(1)));
df_dD2 = 1/(4*pi*sqrt(D(2)));

for i=1:17       % beam %
 dg7_da(i) =  1/20 * df_dD1*(V(1:6,i)'*(dk_da_beam(:,:,i) - D(1:6,1:6)*dM_da_beam(:,:,i)))*V(1:6,i);
end
dg7_da

plot(A, dg1_da)





