%=========================================================================
%
%      FEM code for Damage Consitutive Model (one element test)
%
%=========================================================================
function [r_p,r_q,r_eps,r_epsel,r_epsed,r_epsE,r_epsid,r_Omega,r_f,r_sigmaT] = DSID_v2()

clc
close all
clear all

global c h nu0 E0
global a1 a2 a3 a4
global C0 C1 phi alpha cc
global FTOL

% Initial state

% shao's model
E0 = 6.8E10;  % unit=Pa
nu0 = 0.21;
c = -nu0;
h = 16*(1-nu0^2)/3/E0/(2-nu0);


a1 = -c*h/70;
a2 = (7+2*c)*h/7;
a3 = c*h/7;
a4 = -c*h/35;

% Hayakawa's model

%E0 = 169E9;  % unit=Pa
%nu0 = 0.285;
%a1 = -3.95e-16;
%a2 = 2.5e-12;
%a3 = -4e-13;
%a4 = 4e-12;

phi = 30/180*pi; 
cc=110000;
alpha =2*sin(phi)/sqrt(3)/(3-sin(phi));

C0 = 1.1E5;
%C0 = 273000; %Hayakawa's
C1 = 2.2E6;

Pars=zeros(1,9);
Pars(1)=a1;
Pars(2)=a2;
Pars(3)=a3;
Pars(4)=a4;
Pars(5)=C0;
Pars(6)=C1;
Pars(7)=alpha;
Pars(8)=E0;
Pars(9)=nu0;

FTOL = 1e-6;
Iter=150;

% Load path:
%  1 triaxial compression
% eta = [0 3 3]; % 0 = iso, 3 = triax , 2 = pure shear, 4 = normal and
% shear stress

% Loading:
eta = [0 3];
steps = 2; %steps
chargEps1 = [0 0.1];      
chargs    = [0 0]; % 
ninc = [10 1000];


% %  2 pure shear
% eta = [0 2]; % 0 = iso, 3 = triax , 2 = pure shear
% steps = 2; %steps
% chargEps1 = [0 0];      
% chargs    = [0 100e6]; % 
% ninc = [0 1000];
% %  3 cyclic loading
% eta = [0 3 3 3 3]; % 0 = iso, 3 = triax , 2 = pure shear


% % loading
% eta = [0 4 4];
% steps = 3;
% chargEps1 = [0 0 0];
% chargs = [0 100e6 0;   % s11 vertical
%           0 10e6 0;   % s22
%           0 10e6 0;    % s33
%           0 0 60e6;    % s12
%           0 0 0;    % s23
%           0 0 0];   % s13
% ninc = [10 1000 1000];



% % Loading:
% steps = 5; %steps
% chargEps1 = [0 -0.0005 0.002 -0.005 0.008];      
% chargs    = [0 0 0 0 0]; % 
% ninc = [0 1000 1000 1000 1000];
%
%========================================================================
%         INITIALIZATIONS:
%========================================================================
sigmaT = zeros(3,3); % non-load
sigmaT_v = mat2_mat1(sigmaT);
Omega = zeros(3,3); % undamaged
Omega_v = mat2_mat1(Omega);  
[De0,S0] = matDO1(Omega);   % undamaded stiffness

%========================================================================
%         Strain decomposition:
%         epsT = epsel + epsed + epsid
%========================================================================
epsT = zeros(3,3);
epsel = zeros(3,3); 
epsed = zeros(3,3);
epsE = zeros(3,3);
epsid = zeros(3,3);
% transfer matrix to vector
epsT_v = mat2_mat1(epsT); 
epsel_v = mat2_mat1(epsel); 
epsed_v = mat2_mat1(epsed);
epsE_v = mat2_mat1(epsE);
epsid_v = mat2_mat1(epsid);

% Storage
r_p(1) = trace(sigmaT)/3;%+pin;
r_q(1) = sigmaT(1,1)-sigmaT(2,2);
r_sigmaT(1,:) = sigmaT_v;
r_epsel(1,:) = epsel_v;
r_epsed(1,:) = epsed_v;
r_epsid(1,:) = epsid_v;
r_epsE(1,:) = epsE_v;
r_Omega(1,:) = Omega_v;

% yield function
r_f(1) = fdDP(sigmaT,zeros(3,3),Omega,Pars); %trial test

if r_f(1) > 0
    error('non-elastic initial state')
end

%========================================================================
% SIMULATION (for a single stress path component):
%========================================================================

tinc = 1;
for icharg = 1:steps 
    disp(['============= load step #',num2str(icharg),' ============='])
    
    for inc = 1:ninc(icharg) % load increments
        disp(['              increments #',num2str(inc),'              '])
        
        [matDOm,SOm] = matDO1(Omega);% matDOm SOm 4th tensor
        
        if eta(icharg) == 0 %iso
            if chargs(icharg)~=0   % stress controlled
                dsig = chargs(icharg)/ninc(icharg)*eye(3);
                deps = Aijkl_Bkl(SOm,dsig);
            else    % strain controlled
                deps = chargEps1(icharg)/ninc(icharg)*eye(3);
                dsig = Aijkl_Bkl(matDOm,deps);
                
            end
            
        elseif eta(icharg) == 3 %triaxial
            % Elastic trial:
            if chargs(icharg)~=0   % stress controlled
                dsig = chargs(icharg)/ninc(icharg)*[1 0 0;0 0 0;0 0 0];
                deps = Aijkl_Bkl(SOm,dsig);
            else    % strain controlled
                deps = chargEps1(icharg)/ninc(icharg)*[1 0 0;0 0 0;0 0 0];
                matDOm_2=mat4_mat2(matDOm);% matDOm_2 matrix
                
                deps(2,2) = -matDOm_2(2,1)*deps(1,1)/(matDOm_2(2,2)+matDOm_2(2,3));
                deps(3,3) = deps(2,2);
                dsig = Aijkl_Bkl(matDOm,deps);
                dsig(2,2)=0;
                dsig(3,3)=0;
             
            end
        elseif eta(icharg) == 2 % shear
            % Elastic trial:
            if chargs(icharg)~=0   % stress controlled
                dsig = chargs(icharg)/ninc(icharg)*[0 1 0;1 0 0;0 0 0];
                deps = Aijkl_Bkl(SOm,dsig);
            else    % strain controlled
                deps = chargEps1(icharg)/ninc(icharg)*[0 1 0;1 0 0;0 0 0];
                matDOm_2=mat4_mat2(matDOm);% matDOm_2 matrix
                
                deps(2,2) = -matDOm_2(2,1)*deps(1,1)/(matDOm_2(2,2)+matDOm_2(2,3));
                deps(3,3) = deps(2,2);
                dsig = Aijkl_Bkl(matDOm,deps);
                dsig(2,2)=0;
                dsig(3,3)=0;
             
            end
        elseif eta(icharg) == 4 % shear and normal stress
             % only stress controlled
             dsig = 1/ninc(icharg)*[chargs(1,icharg) chargs(4,icharg) chargs(6,icharg);
                                    chargs(4,icharg) chargs(2,icharg) chargs(5,icharg);
                                    chargs(6,icharg) chargs(5,icharg) chargs(3,icharg)];
             deps = Aijkl_Bkl(SOm,dsig);
        end
               
        fd0 = fdDP(sigmaT,zeros(3,3),Omega,Pars);
        fd = fdDP(sigmaT,dsig,Omega,Pars);
% 
        nbStep = 0;
        if fd <= FTOL % elastic increment

            sigmaT = sigmaT + dsig;
            %deps = Aijkl_Bkl(SOm,dsig);
            depsel = Aijkl_Bkl(S0,dsig);
            depsid = zeros(3,3);
            depsE = deps;
            depsed = deps - depsel;

            epsE = epsE + depsE;
            epsT = epsT + deps ;
            epsed = epsed + depsed;
            epsid = epsid + depsid;
            epsel = epsel + depsel; 

            pn = trace(sigmaT)/3;
            qn = sigmaT(1,1) - sigmaT(2,2);
            
            end_r = length(r_p)+1;
            r_p(end_r) = pn;
            r_q(end_r) = qn;
            % transfer matrix to vector
            sigmaT_v = mat2_mat1(sigmaT);
            epsT_v = mat2_mat1(epsT);
            epsel_v = mat2_mat1(epsel); 
            epsed_v = mat2_mat1(epsed);
            epsE_v = mat2_mat1(epsE);
            epsid_v = mat2_mat1(epsid);
            Omega_v = mat2_mat1(Omega);
            r_sigmaT(end_r,:)=sigmaT_v;
            r_eps(end_r,:) = epsT_v;
            r_epsel(end_r,:) = epsel_v;
            r_epsed(end_r,:) = epsed_v;
            r_epsE(end_r,:) = epsE_v;
            r_epsid(end_r,:) = epsid_v;
            r_Omega(end_r,:) = Omega_v;
            r_f(end_r) = fdDP(sigmaT,zeros(3,3),Omega,Pars);
            
        else
           
            depsid = zeros(3,3);
            depsE = zeros(3,3);
            if (chargs(icharg)~=0 || eta(icharg) == 4)
                 incinc=1;
                lambdam=0;lambdan=10^-3;
              while (abs(fd)>FTOL) && (incinc<Iter) %  for incinc=1:50%
               [ fdm,epsm,Omegam,depsidm ] = fd_lam2( lambdam, Omega, sigmaT, epsid, Pars,dsig);
               [ fdn,epsn,Omegan,depsidn ] = fd_lam2( lambdan, Omega, sigmaT, epsid, Pars,dsig);
               lambda = lambdan - fdn*(lambdan-lambdam)/(fdn-fdm);
                [fd,eps,OmegaT,depsidT] = fd_lam2( lambda, Omega, sigmaT, epsid, Pars,dsig);          
                dOmega = OmegaT-Omega;                
                sigmaTT = sigmaT + dsig;
                
                fd = fdDP(sigmaTT,zeros(3,3),OmegaT,Pars);
                incinc=incinc+1;
                lambdam=lambdan;
                lambdan=lambda;
              end
                deps = eps - epsT;
                depsid = depsidT;
                depsel = Aijkl_Bkl(S0,dsig);
                depsed = deps-depsid-depsel;
                depsE = depsed + depsel;
                epsT = epsT + deps;
                Omega = Omega + dOmega;

%

                depsed = depsE - depsel;
                epsel = epsel + depsel; 
                epsid = epsid + depsid;
                epsE = epsE + depsE;  
                epsed = epsed + depsed;
%                 fd = fdDP(sigmaTT,zeros(3,3),Omega,Pars);
%                end
%                pause

            else
                incinc=1;
                lambdam=0;lambdan=10^-3;
              while (abs(fd)>FTOL) && (incinc<Iter) %  for incinc=1:50%
                [fdm,sigm,Omegam,depsidm] = fd_lam( lambdam, Omega, sigmaT, epsT, epsid, Pars,deps);
                [fdn,sign,Omegan,depsidn] = fd_lam( lambdan, Omega, sigmaT, epsT, epsid, Pars,deps);
                lambda = lambdan - fdn*(lambdan-lambdam)/(fdn-fdm);
                [fd,sig,OmegaT,depsidT] = fd_lam( lambda, Omega, sigmaT, epsT, epsid, Pars,deps);
                sigmaTT = sig;          
                dsig = sigmaTT-sigmaT;
                dOmega = OmegaT-Omega;
                fd = fdDP(sigmaTT,zeros(3,3),OmegaT,Pars);
                incinc=incinc+1;
                lambdam=lambdan;
                lambdan=lambda;
              end
              %pause
              dY_dsig = dY_dsigf(sigmaT+dsig);
              depsid=depsidT;
              epsT = epsT + deps;
              Omega = Omega + dOmega;
              depsel = Aijkl_Bkl(S0,dsig);
              depsed = deps - depsel-depsid;
              depsE = depsel+depsed;
              
              epsel = epsel + depsel; 
              epsid = epsid + depsid;
              epsE = epsE + depsE;       
              epsed = epsed + depsed;           
            end
            sigmaT = sigmaT + dsig;
            pn = trace(sigmaT)/3;
            qn = sigmaT(1,1) - sigmaT(2,2);
            end_r = length(r_p)+1;
            r_p(end_r) = pn;
            r_q(end_r) = qn;
            sigmaT_v=mat2_mat1(sigmaT);
            epsT_v = mat2_mat1(epsT);
            epsel_v = mat2_mat1(epsel); 
            epsed_v = mat2_mat1(epsed);
            epsE_v = mat2_mat1(epsE);
            epsid_v = mat2_mat1(epsid);
            Omega_v = mat2_mat1(Omega);
            r_sigmaT(end_r,:)=sigmaT_v;
            r_eps(end_r,:) = epsT_v;
            r_epsel(end_r,:) = epsel_v;
            r_epsed(end_r,:) = epsed_v;
            r_epsE(end_r,:) = epsE_v;
            r_epsid(end_r,:) = epsid_v;
            r_Omega(end_r,:) = Omega_v; 

%         fd
         fd1=fdDP(sigmaT,zeros(3,3),Omega,Pars);
            r_f(end_r) = fdDP(sigmaT,zeros(3,3),Omega,Pars);
            


        end
        tinc = tinc + 1;% just counts total increments (if several loadings)
    end % increments
    r_ends(icharg) = length(r_p);
     
end % loading parts
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST-PROCESSING:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

figure('Name','q(eps1)','NumberTitle','off')
plot(r_eps(:,1),r_q/1000000,'Linewidth',3)
xlabel('Axial strain, \epsilon_1','FontSize',20)
ylabel('Deviatoric stress, q (MPa)','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'q_eps1.eps'

figure1=figure('Name','q(eps1,eps3)','NumberTitle','off');
plot(r_eps(:,1),r_q/1000000,'-b',r_eps(:,3),r_q/1000000,'-r','Linewidth',3)
xlabel('Axial strain, \epsilon_1','FontSize',20)
ylabel('Deviatoric stress, q (MPa)','FontSize',20)
legend('\epsilon_1', '\epsilon_3','Location','Best')
grid
set(gca,'FontSize',20)
set(gca,'XDir','reverse')
set(gca,'YDir','reverse')
print -depsc 'q_eps1-3.eps'

figure('Name','q(eps3)','NumberTitle','off')
plot(r_eps(:,3),r_q/1000000,'Linewidth',3)
xlabel('Axial strain, \epsilon_3','FontSize',20)
ylabel('Deviatoric stress, q (MPa)','FontSize',20)

grid
set(gca,'FontSize',20)
print -depsc 'q_eps3.eps'

figure('Name','Omega(q)','NumberTitle','off')
plot(r_q/1000000,r_Omega(:,1),'-b','Linewidth',3)
hold on
plot(r_q/1000000,r_Omega(:,3),'--r','Linewidth',5)
xlabel('Deviatoric stress, q (MPa)','FontSize',20)
ylabel('Damage variable, D','FontSize',20)
legend('D_{11}', 'D_{33}','Location','Best')
grid
set(gca,'FontSize',20)
print -depsc 'Omega_q.eps'


%&&&&&&&&&

figure('Name','Omega(eps1)','NumberTitle','off')
plot(r_eps(:,1),r_Omega(:,1),'-b','Linewidth',3)
hold on
plot(r_eps(:,1),r_Omega(:,3),'--r','Linewidth',5)
xlabel('Axial strain, \epsilon_1','FontSize',20)
ylabel('Damage variable, D','FontSize',20)
legend('D_{11}', 'D_{33}','Location','Best')
grid
set(gca,'FontSize',20)
print -depsc 'Omega_eps1.eps'


figure('Name','sigma12(eps12)','NumberTitle','off')
plot(r_eps(:,4),r_sigmaT(:,4)/1000000,'Linewidth',3)
xlabel('Axial strain, \epsilon_{12}','FontSize',20)
ylabel('Deviatoric stress, \sigma_{12} (MPa)','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'sigma12_eps12.eps'

figure('Name','Omega3(eps3)','NumberTitle','off')
plot(r_eps(:,3),r_Omega(:,3),'Linewidth',3)
xlabel('Axial strain, \epsilon_3','FontSize',20)
ylabel('Damage variable, \Omega_3 (-)','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'Omega3_eps3.eps'

figure('Name','Omega1(eps1)','NumberTitle','off')
plot(r_eps(:,1),r_Omega(:,1),'Linewidth',3)
xlabel('Axial strain, \epsilon_1','FontSize',20)
ylabel('Damage variable, \Omega_1 (-)','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'Omega1_eps1.eps'

figure('Name','Omega3(eps1)','NumberTitle','off')
plot(r_eps(:,1),r_Omega(:,3),'Linewidth',3)
xlabel('Axial strain, \epsilon_1','FontSize',20)
ylabel('Damage variable, \Omega_3 (-)','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'Omega3_eps1.eps'

figure('Name','Omega2(eps1)','NumberTitle','off')
plot(r_eps(:,1),r_Omega(:,2),'Linewidth',3)
xlabel('Axial strain, \epsilon_1','FontSize',20)
ylabel('Damage variable, \Omega_2 (-)','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'Omega2_eps1.eps'

figure('Name','Omega1(eps12)','NumberTitle','off')
plot(r_eps(:,4),r_Omega(:,1),'Linewidth',3)
xlabel('Axial strain, \epsilon_{12}','FontSize',20)
ylabel('Damage variable, \Omega_1 (-)','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'Omega1_eps12.eps'

figure('Name','Omega12(eps12)','NumberTitle','off')
plot(r_eps(:,4),r_Omega(:,4),'Linewidth',3)
xlabel('Axial strain, \epsilon_{12}','FontSize',20)
ylabel('Damage variable, \Omega_{12} (-)','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'Omega12_eps12.eps'

figure('Name','f','NumberTitle','off')
plot(r_f,'Linewidth',3)
xlabel('number','FontSize',20)
ylabel('Damage function, f','FontSize',20)
grid
set(gca,'FontSize',20)
print -depsc 'f.eps'


end



%==========================================================================
%
%                  Functions
%
%==========================================================================

function [ scalar ] = Aij_Bij( A,B )
%   Double contraction: scalar = A_ij*B_ij

scalar = 0;
for i= 1:3
    for j=1:3 
        scalar = scalar + A(i,j)*B(i,j); 
    end
end
   
end

function [ C ] = Aij_Bkl( A,B )

C(1:3,1:3,1:3,1:3) = 0;
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                C(i,j,k,l) = A(i,j)*B(k,l);
            end
        end
    end
end

end

function [ C ] = Aijkl_Bij( A,B )

C = zeros(3,3);
for k = 1:3
    for l = 1:3
        for i = 1:3
            for j = 1:3
              C(k,l) = C(k,l)+A(i,j,k,l)*B(i,j);
            end
        end
    end
end

end

function [ C ] = Aijkl_Bkl( A,B )

C = zeros(3,3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
              C(i,j) = C(i,j)+A(i,j,k,l)*B(k,l);
            end
        end
    end
end

end

function [ depsid ] = depsidf(sigmaT,dsig,invA)
%   Calculate irreversible strain increment
global a1 a2 a3 a4
global C0 C1 phi alpha cc

E = eye(3);
temp1=zeros(3,3);
temp2=zeros(3,3);
temp3=zeros(3,3);
dY_dsig = dY_dsigf(sigmaT);
dY = Aijkl_Bkl(dY_dsig,dsig);
Yd1 = a1*(trace(sigmaT))^2*E+a2*sigmaT*sigmaT+a3*trace(sigmaT)*sigmaT+a4*trace(sigmaT*sigmaT)*E;
P_1 = matP_1(sigmaT);
P_2 = matP_2(sigmaT);
F1ij = Aijkl_Bkl(P_1,Yd1)-1/3*Aij_Bij(Aijkl_Bkl(P_1,Yd1),E)*E;
F2ij = Aijkl_Bkl(P_2,Yd1);
P_3 = P_1-1/3*Aij_Bkl(E,Aijkl_Bij(P_1,E));
df_dOmega = -C1*E;
dg_dY = Aijkl_Bij(P_2,F2ij)/sqrt(2*Aij_Bij(F2ij,F2ij));
df_dY = Aijkl_Bij(P_3,F1ij)/sqrt(2*Aij_Bij(F1ij,F1ij))-alpha*Aijkl_Bij(P_1,E);
df_dsig = Aijkl_Bij(dY_dsig,df_dY);
temp1= Aijkl_Bkl(invA,df_dY);% A^(-1):df_dY
temp2=Aijkl_Bij(dY_dsig,temp1);% A^(-1):df_dY:dY_dsig
temp3=Aijkl_Bij(invA,df_dsig);% df_dsig:A^(-1)
temp6=Aij_Bij(temp3,dsig);%df_dsig:A^(-1):dsig
temp4=Aij_Bij(df_dsig,temp2);%df_dsig:A^(-1):df_dY:dY_dsig
temp5=Aij_Bij(df_dOmega,dg_dY);%df_dOmega:dg_dY

lambdad=temp6/(temp4-temp5);
depsid = lambdad*Aijkl_Bij(dY_dsig,df_dY);  %% dg_dY ===> df_dY
end

function [ dY_dsig ] = dY_dsigf( sigmaT )
%
global E0 nu0 a1 a2 a3 a4
E = eye(3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                dY_dsig(i,j,k,l)=2*a1*trace(sigmaT)*E(i,j)*E(k,l)+1/2*a2*(E(i,k)*sigmaT(l,j)+...
                    E(i,l)*sigmaT(j,k)+E(j,l)*sigmaT(i,k)+E(j,k)*sigmaT(i,l))+a3*(E(k,l)*sigmaT(i,j)+...
                1/2*trace(sigmaT)*(E(i,k)*E(j,l)+E(i,l)*E(j,k)))+2*a4*sigmaT(k,l)*E(i,j);
            end
        end
    end
end

end

function [ f ] = fdDP( sigmaT,dsig,Omega,Pars )
%=========================================================================
%
%       Damage Yield Function (Pseudo-Drucker-Pager)
%
%=========================================================================
a1=Pars(1);
a2=Pars(2);
a3=Pars(3);
a4=Pars(4);
C0=Pars(5);
C1=Pars(6);
alpha=Pars(7);

sigmaTT = sigmaT + dsig ;
trsigmaT = trace(sigmaTT);
trOmega = trace(Omega);
Yd1 = a1*(trsigmaT)^2*eye(3)+a2*sigmaTT*sigmaTT+a3*trsigmaT*sigmaTT+a4*trace(sigmaTT*sigmaTT)*eye(3);

P_1=matP_1(sigmaTT);

P_1Y=Aijkl_Bkl(P_1,Yd1);
trY = trace(P_1Y);
S = P_1Y-1/3*trY*eye(3);
f = sqrt(0.5*Aij_Bij(S,S))-alpha*trY-C0-C1*trOmega;
end

function [ invA ] = invmat4( A )
% calculate inverse of 4th-order tensor
E6=eye(6);
A_2=mat4_mat2(A);
invA_2=A_2\E6;
invA=mat2_mat4(invA_2);

end

function [ vector ] = mat2_mat1( matrix )
%==========================================================================
%
%    MAT2_MAT1 
%    Transfer a 3*3 matrix to 6*1 vextor (default format in ABAQUS)
%
%    sig_11 sig_12 sig_13
%    sig_21 sig_22 sig_23 ==> [sig_11 sig_22 sig_33 sig_12 sig_13 sig_23]'
%    sig_31 sig_32 sig_33
%
%==========================================================================
vector = zeros(6,1);
for i = 1:3
    vector(i) = matrix(i,i);
    for j = i+1:3
        vector(i+j+1) = matrix(i,j);
    end
end

end

function [ tensor ] = mat2_mat4( matrix )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
      tensor(1:3,1:3,1:3,1:3) = 0;
      
      tensor(1,1,1,1) = matrix(1,1);
      tensor(1,1,2,2) = matrix(1,2);
      tensor(1,1,3,3) = matrix(1,3);
      tensor(1,1,1,2) = matrix(1,4);
      tensor(1,1,2,1) = matrix(1,4);
      tensor(1,1,2,3) = matrix(1,5);
      tensor(1,1,3,2) = matrix(1,5);
      tensor(1,1,1,3) = matrix(1,6);
      tensor(1,1,3,1) = matrix(1,6);

      tensor(2,2,1,1) = matrix(2,1);
      tensor(2,2,2,2) = matrix(2,2);
      tensor(2,2,3,3) = matrix(2,3);
      tensor(2,2,1,2) = matrix(2,4);
      tensor(2,2,2,1) = matrix(2,4);
      tensor(2,2,2,3) = matrix(2,5);
      tensor(2,2,3,2) = matrix(2,5);
      tensor(2,2,1,3) = matrix(2,6);
      tensor(2,2,3,1) = matrix(2,6);

      tensor(3,3,1,1) = matrix(3,1);
      tensor(3,3,2,2) = matrix(3,2);
      tensor(3,3,3,3) = matrix(3,3);
      tensor(3,3,1,2) = matrix(3,4);
      tensor(3,3,2,1) = matrix(3,4);
      tensor(3,3,2,3) = matrix(3,5);
      tensor(3,3,3,2) = matrix(3,5);
      tensor(3,3,1,3) = matrix(3,6);
      tensor(3,3,3,1) = matrix(3,6);

      tensor(1,2,1,1) = matrix(4,1);
      tensor(1,2,2,2) = matrix(4,2);
      tensor(1,2,3,3) = matrix(4,3);
      tensor(1,2,1,2) = matrix(4,4);
      tensor(1,2,2,1) = matrix(4,4);
      tensor(1,2,2,3) = matrix(4,5);
      tensor(1,2,3,2) = matrix(4,5);
      tensor(1,2,1,3) = matrix(4,6);
      tensor(1,2,3,1) = matrix(4,6);

      tensor(2,3,1,1) = matrix(5,1);
      tensor(2,3,2,2) = matrix(5,2);
      tensor(2,3,3,3) = matrix(5,3);
      tensor(2,3,1,2) = matrix(5,4);
      tensor(2,3,2,1) = matrix(5,4);
      tensor(2,3,2,3) = matrix(5,5);
      tensor(2,3,3,2) = matrix(5,5);
      tensor(2,3,1,3) = matrix(5,6);
      tensor(2,3,3,1) = matrix(5,6);

      tensor(1,3,1,1) = matrix(6,1);
      tensor(1,3,2,2) = matrix(6,2);
      tensor(1,3,3,3) = matrix(6,3);
      tensor(1,3,1,2) = matrix(6,4);
      tensor(1,3,2,1) = matrix(6,4);
      tensor(1,3,2,3) = matrix(6,5);
      tensor(1,3,3,2) = matrix(6,5);
      tensor(1,3,1,3) = matrix(6,6);
      tensor(1,3,3,1) = matrix(6,6);
      
      tensor(2,1,1,1) = matrix(4,1);
      tensor(2,1,2,2) = matrix(4,2);
      tensor(2,1,3,3) = matrix(4,3);
      tensor(2,1,1,2) = matrix(4,4);
      tensor(2,1,2,1) = matrix(4,4);
      tensor(2,1,2,3) = matrix(4,5);
      tensor(2,1,3,2) = matrix(4,5);
      tensor(2,1,1,3) = matrix(4,6);
      tensor(2,1,3,1) = matrix(4,6);

      tensor(3,2,1,1) = matrix(5,1);
      tensor(3,2,2,2) = matrix(5,2);
      tensor(3,2,3,3) = matrix(5,3);
      tensor(3,2,1,2) = matrix(5,4);
      tensor(3,2,2,1) = matrix(5,4);
      tensor(3,2,2,3) = matrix(5,5);
      tensor(3,2,3,2) = matrix(5,5);
      tensor(3,2,1,3) = matrix(5,6);
      tensor(3,2,3,1) = matrix(5,6);

      tensor(3,1,1,1) = matrix(6,1);
      tensor(3,1,2,2) = matrix(6,2);
      tensor(3,1,3,3) = matrix(6,3);
      tensor(3,1,1,2) = matrix(6,4);
      tensor(3,1,2,1) = matrix(6,4);
      tensor(3,1,2,3) = matrix(6,5);
      tensor(3,1,3,2) = matrix(6,5);
      tensor(3,1,1,3) = matrix(6,6);
      tensor(3,1,3,1) = matrix(6,6);

      
      
end

function [ matrix ] = mat4_mat2( tensor )

      matrix = zeros (6,6);

      matrix(1,1)=tensor(1,1,1,1);
      matrix(1,2)=tensor(1,1,2,2);
      matrix(1,3)=tensor(1,1,3,3);
      matrix(1,4)=tensor(1,1,1,2);
      matrix(1,5)=tensor(1,1,2,3);
      matrix(1,6)=tensor(1,1,1,3);

      matrix(2,1)=tensor(2,2,1,1);
      matrix(2,2)=tensor(2,2,2,2);
      matrix(2,3)=tensor(2,2,3,3);
      matrix(2,4)=tensor(2,2,1,2);
      matrix(2,5)=tensor(2,2,2,3);
      matrix(2,6)=tensor(2,2,1,3);

      matrix(3,1)=tensor(3,3,1,1);
      matrix(3,2)=tensor(3,3,2,2);
      matrix(3,3)=tensor(3,3,3,3);
      matrix(3,4)=tensor(3,3,1,2);
      matrix(3,5)=tensor(3,3,2,3);
      matrix(3,6)=tensor(3,3,1,3);

      matrix(4,1)=tensor(1,2,1,1);
      matrix(4,2)=tensor(1,2,2,2);
      matrix(4,3)=tensor(1,2,3,3);
      matrix(4,4)=tensor(1,2,1,2);
      matrix(4,5)=tensor(1,2,2,3);
      matrix(4,6)=tensor(1,2,1,3);

      matrix(5,1)=tensor(2,3,1,1);
      matrix(5,2)=tensor(2,3,2,2);
      matrix(5,3)=tensor(2,3,3,3);
      matrix(5,4)=tensor(2,3,1,2);
      matrix(5,5)=tensor(2,3,2,3);
      matrix(5,6)=tensor(2,3,1,3);

      matrix(6,1)=tensor(1,3,1,1);
      matrix(6,2)=tensor(1,3,2,2);
      matrix(6,3)=tensor(1,3,3,3);
      matrix(6,4)=tensor(1,3,1,2);
      matrix(6,5)=tensor(1,3,2,3);
      matrix(6,6)=tensor(1,3,1,3);


end

function [C] = Aijpq_Bpqkl(A,B)

C = zeros(3,3,3,3);

  for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for p=1:3
                    for q=1:3
                        C(i,j,k,l)=C(i,j,k,l)+A(i,j,p,q)*B(p,q,k,l);
                    end
                end
            end
        end
    end
  end

end

function [ C ] = Aijklpq_Bpq(A,B)

C = zeros(3,3,3,3);

  for i=1:3
      for j=1:3
        for k=1:3
            for l=1:3
                for p=1:3
                    for q=1:3
                        C(i,j,k,l)=C(i,j,k,l)+A(i,j,k,l,p,q)*B(p,q);
                    end
                end
            end
        end
    end
  end

end

function [ matinvA ] = matA(sigmaT,SOm1,H)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

global E0 nu0 a1 a2 a3 a4
b1 = (1+nu0)/E0/2;
b2 = nu0/E0;
E = eye(3);
%E6 = eye(6);
%trOmega = trace(Omega);
dS_dOmega(1:3,1:3,1:3,1:3,1:3,1:3) =0;
sig_dsdOmega(1:3,1:3,1:3,1:3)=0;
temp1(1:3,1:3,1:3,1:3)=0;
temp2(1:3,1:3,1:3,1:3)=0;
matinvA(1:3,1:3,1:3,1:3)=0;
dY_dsig=dY_dsigf(sigmaT);

for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for p=1:3
                    for q=1:3
                        dS_dOmega(i,j,k,l,p,q)=2*a1*E(i,j)*E(k,l)*E(p,q)+...
                            1/4*a2*(E(i,k)*(E(p,j)*E(q,l)+E(p,l)*E(q,j))+...
                        E(i,l)*(E(p,j)*E(q,k)+E(p,k)*E(q,j))+...
                        E(j,l)*(E(i,p)*E(q,k)+E(i,q)*E(p,k))+...
                        E(j,k)*(E(i,p)*E(q,l)+E(i,q)*E(p,l)))+...
                        1/2*a3*(E(i,j)*(E(k,p)*E(l,q)+E(k,q)*E(l,p))+E(k,l)*(E(i,p)*E(j,q)+E(i,q)*E(j,p)))+...
                        a4*(E(i,k)*E(j,l)+E(i,l)*E(j,k))*E(p,q);
                    end
                end
            end
        end
    end
end


for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for p=1:3
                    for q=1:3
                        sig_dsdOmega(i,j,k,l)=sig_dsdOmega(i,j,k,l)+sigmaT(p,q)*dS_dOmega(p,q,i,j,k,l);
                    end
                end
            end
        end
    end
end


temp1 = Aijpq_Bpqkl(sig_dsdOmega,H);
temp2 = Aijpq_Bpqkl(temp1,dY_dsig);

matinvA = temp2+SOm1;
end

function [ matCed ] = matCedf( sigmaT,invA )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global a1 a2 a3 a4 C0 C1 alpha

E = eye(3);
temp1=zeros(3,3);
temp2=zeros(3,3);
temp3=zeros(3,3);
Yd1 = a1*(trace(sigmaT))^2*E+a2*sigmaT*sigmaT+a3*trace(sigmaT)*sigmaT+a4*trace(sigmaT*sigmaT)*E;
dY_dsig = dY_dsigf(sigmaT);
P_1 = matP_1(sigmaT);
P_2 = matP_2(sigmaT);
F1ij = Aijkl_Bkl(P_1,Yd1)-1/3*Aij_Bij(Aijkl_Bkl(P_1,Yd1),E)*E;
F2ij = Aijkl_Bkl(P_2,Yd1);
df_dOmega = -C1*E;
dg_dY = Aijkl_Bij(P_2,F2ij)/sqrt(2*Aij_Bij(F2ij,F2ij));
P_3 = P_1-1/3*Aij_Bkl(E,Aijkl_Bij(P_1,E));
df_dY = Aijkl_Bij(P_3,F1ij)/sqrt(2*Aij_Bij(F1ij,F1ij))-alpha*Aijkl_Bij(P_1,E);
%pause
df_dsig = Aijkl_Bij(dY_dsig,df_dY);
%matH=-Aij_Bkl(dg_dY,df_dY)/Aij_Bij(df_dOmega,dg_dY);
temp1= Aijkl_Bkl(invA,df_dY); %===>A^(-1):df_dY
temp2=Aijkl_Bij(dY_dsig,temp1);%A^(-1):dg_dY:dY_dsig
temp3=Aijkl_Bij(invA,df_dsig);%df_dsig:A^(-1)
temp4=Aij_Bij(df_dsig,temp2);%df_dsig:A^(-1):dg_dY:dY_dsig
temp5=Aij_Bij(df_dOmega,dg_dY);%df_dOmega:dg_dY
matCed=invA-Aij_Bkl(temp2,temp3)/(temp4-temp5);
end

function [matDz, matS] = matDO1(Omega)
global E0 nu0 a1 a2 a3 a4
b1 = (1+nu0)/E0/2;
b2 = nu0/E0;
E = eye(3);
E6 = eye(6);
trOmega = trace(Omega);
matS(1:3,1:3,1:3,1:3) = 0;
matS_2(1:6,1:6) = 0;
matDz(1:3,1:3,1:3,1:3) = 0;
matDz_2(1:6,1:6) = 0;
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                matS(i,j,k,l) = b1*(E(i,k)*E(j,l)+E(i,l)*E(j,k))-...
                    b2*E(i,j)*E(k,l)+2*a1*trOmega*E(i,j)*E(k,l)+...
                    0.5*a2*(E(i,k)*Omega(j,l)+E(i,l)*Omega(j,k)+...
                    Omega(i,k)*E(j,l)+Omega(i,l)*E(j,k))+...
                    a3*(E(i,j)*Omega(k,l)+Omega(i,j)*E(k,l))+...
                    a4*trOmega*(E(i,k)*E(j,l)+E(i,l)*E(j,k));
            end
        end
    end
end    
matS_2 = mat4_mat2(matS);
matDz_2 = matS_2\E6;
matDz = mat2_mat4(matDz_2);


end

function [ matH ] = matH(sigmaT)
%UNTITLED3 Summary of this function goes here
%   matH_ijkl
%   dOmega = matH_ijkl*dY_kl
global a1 a2 a3 a4
global C0 C1 phi alpha cc

matH(1:3,1:3,1:3,1:3)=0;
E = eye(3);
Yd1 = a1*(trace(sigmaT))^2*E+a2*sigmaT*sigmaT+a3*trace(sigmaT)*sigmaT+a4*trace(sigmaT*sigmaT)*E;
P_1 = matP_1(sigmaT);
P_2 = matP_2(sigmaT);

F1ij = Aijkl_Bkl(P_1,Yd1)-1/3*Aij_Bij(Aijkl_Bkl(P_1,Yd1),E)*E;
F2ij = Aijkl_Bkl(P_2,Yd1);

df_dOmega = -C1*E;
dg_dY = Aijkl_Bij(P_2,F2ij)/sqrt(2*Aij_Bij(F2ij,F2ij));
P_3 = P_1-1/3*Aij_Bkl(E,Aijkl_Bij(P_1,E));
df_dY = Aijkl_Bij(P_3,F1ij)/sqrt(2*Aij_Bij(F1ij,F1ij))-alpha*Aijkl_Bij(P_1,E);
matH=-Aij_Bkl(dg_dY,df_dY)/Aij_Bij(df_dOmega,dg_dY);



end

function [ P_1 ] = matP_1(sigmaT)
%UNTITLED4 Summary of this function goes here
%   P_1 projection tensor
%
%   P_1=H(sigmaT)-H(-sigmaT)
%
%
P_1(1:3,1:3,1:3,1:3)=0;

n1=[1;0;0];
n2=[0;1;0];
n3=[0;0;1];
[V,D] = eig(sigmaT);
n1=V(:,1);
n2=V(:,2);
n3=V(:,3);
P_1=(heaviside(D(1,1))-heaviside(-D(1,1)))*ni_nj_nk_nl(n1)+...
    (heaviside(D(2,2))-heaviside(-D(2,2)))*ni_nj_nk_nl(n2)+...
    (heaviside(D(3,3))-heaviside(-D(3,3)))*ni_nj_nk_nl(n3);

end

function [ P_2 ] = matP_2(sigmaT)
%UNTITLED4 Summary of this function goes here
%   P_1 projection tensor
%
%   P_1=H(sigmaT)-H(-sigmaT)
%
%
P_2(1:3,1:3,1:3,1:3)=0;

n1=[1;0;0];
n2=[0;1;0];
n3=[0;0;1];
[V,D] = eig(sigmaT);
n1=V(:,1);
n2=V(:,2);
n3=V(:,3);
s(1:3)=0;
s(1)=D(1,1)-max(D(1,1),max(D(2,2),D(3,3)));
s(2)=D(2,2)-max(D(1,1),max(D(2,2),D(3,3)));
s(3)=D(3,3)-max(D(1,1),max(D(2,2),D(3,3)));
for ii=1:3
    if s(ii)==0
        s(ii)=0;
    else
        s(ii)=heaviside(-s(ii));
    end
end

P_2=s(1)*ni_nj_nk_nl(n1)+...
    s(2)*ni_nj_nk_nl(n2)+...
    s(3)*ni_nj_nk_nl(n3);

end

function [ M ] = ni_nj_nk_nl( n )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
M(1:3,1:3,1:3,1:3)=0;

for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                M(i,j,k,l) = n(i)*n(j)*n(k)*n(l);
            end
        end
    end
end 
end

function [ fdm,sigm,Omegam,depsidm ] = fd_lam( lambdam, Omega, sigmaT, epsT, epsid, Pars,deps)
%    Iteration solving for lambdan
%

a1=Pars(1);
a2=Pars(2);
a3=Pars(3);
a4=Pars(4);
C0=Pars(5);
C1=Pars(6);
alpha=Pars(7);

E=eye(3);
Yd1 = a1*(trace(sigmaT))^2*E+a2*sigmaT*sigmaT+a3*trace(sigmaT)*sigmaT+a4*trace(sigmaT*sigmaT)*E;
dY_dsig = dY_dsigf(sigmaT);
P_1 = matP_1(sigmaT);
P_2 = matP_2(sigmaT);
F1ij = Aijkl_Bkl(P_1,Yd1)-1/3*Aij_Bij(Aijkl_Bkl(P_1,Yd1),E)*E;
F2ij = Aijkl_Bkl(P_2,Yd1);
df_dOmega = -C1*E;
dg_dY = Aijkl_Bij(P_2,F2ij)/sqrt(2*Aij_Bij(F2ij,F2ij));
dY_dsig = dY_dsigf(sigmaT);
P_3 = P_1-1/3*Aij_Bkl(E,Aijkl_Bij(P_1,E));
df_dY = Aijkl_Bij(P_3,F1ij)/sqrt(2*Aij_Bij(F1ij,F1ij))-alpha*Aijkl_Bij(P_1,E);
% updating depsid
depsidm = lambdam*Aijkl_Bij(dY_dsig,df_dY);
% depsidn = lambdan*Aijkl_Bij(dY_dsig,df_dY);
% updating Omega
dOmegam = lambdam*dg_dY;
% dOmegan = lambdan*dg_dY;
Omegam = Omega+dOmegam;
% Omegan = Omega+dOmegan;
[matD,S] = matDO1(Omega);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for p=1:3
                    for q=1:3
                        dS_dOmega(i,j,k,l,p,q)=2*a1*E(i,j)*E(k,l)*E(p,q)+...
                            1/4*a2*(E(i,k)*(E(p,j)*E(q,l)+E(p,l)*E(q,j))+...
                        E(i,l)*(E(p,j)*E(q,k)+E(p,k)*E(q,j))+...
                        E(j,l)*(E(i,p)*E(q,k)+E(i,q)*E(p,k))+...
                        E(j,k)*(E(i,p)*E(q,l)+E(i,q)*E(p,l)))+...
                        1/2*a3*(E(i,j)*(E(k,p)*E(l,q)+E(k,q)*E(l,p))+E(k,l)*(E(i,p)*E(j,q)+E(i,q)*E(j,p)))+...
                        a4*(E(i,k)*E(j,l)+E(i,l)*E(j,k))*E(p,q);
                    end
                end
            end
        end
    end
end
temp1 = -invmat4(Aijpq_Bpqkl(S,S));
dmatD_dOmega=zeros(3,3,3,3,3,3);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for m=1:3
                    for n=1:3
                        for p=1:3
                            for q=1:3
    dmatD_dOmega(i,j,k,l,m,n)=dmatD_dOmega(i,j,k,l,m,n)+temp1(i,j,p,q)*dS_dOmega(p,q,k,l,m,n);
                            end
                        end
                    end
                end
            end
        end
    end
end

matDem = matD + lambdam*Aijklpq_Bpq(dmatD_dOmega,dg_dY);
% matDen = matD + lambdan*Aijklpq_Bpq(dmatD_dOmega,dg_dY);
epsm = epsT+deps-epsid-depsidm;
% epsn = epsT+deps-epsid-depsidn;
sigm = Aijkl_Bkl(matDem,epsm);
% sign = Aijkl_Bkl(matDen,epsn);
% Ydm = a1*(trace(sigm))^2*E+a2*sigm*sigm+a3*trace(sigm)*sigm+a4*trace(sigm*sigm)*E;
% Ydn = a1*(trace(sign))^2*E+a2*sign*sign+a3*trace(sign)*sign+a4*trace(sign*sign)*E;
fdm = fdDP(sigm,zeros(3,3),Omegam,Pars);
% fdn = fdDP(sign,zeros(3,3),Omegan,Pars);
end

function [ fdm,epsm,Omegam,depsidm ] = fd_lam2( lambdam, Omega, sigmaT, epsid, Pars,dsig)
%    Iteration solving for lambdan
%

a1=Pars(1);
a2=Pars(2);
a3=Pars(3);
a4=Pars(4);
C0=Pars(5);
C1=Pars(6);
alpha=Pars(7);

E=eye(3);
Yd1 = a1*(trace(sigmaT))^2*E+a2*sigmaT*sigmaT+a3*trace(sigmaT)*sigmaT+a4*trace(sigmaT*sigmaT)*E;
dY_dsig = dY_dsigf(sigmaT);
P_1 = matP_1(sigmaT);
P_2 = matP_2(sigmaT);
F1ij = Aijkl_Bkl(P_1,Yd1)-1/3*Aij_Bij(Aijkl_Bkl(P_1,Yd1),E)*E;
F2ij = Aijkl_Bkl(P_2,Yd1);
df_dOmega = -C1*E;
dg_dY = Aijkl_Bij(P_2,F2ij)/sqrt(2*Aij_Bij(F2ij,F2ij));
dY_dsig = dY_dsigf(sigmaT);
P_3 = P_1-1/3*Aij_Bkl(E,Aijkl_Bij(P_1,E));
df_dY = Aijkl_Bij(P_3,F1ij)/sqrt(2*Aij_Bij(F1ij,F1ij))-alpha*Aijkl_Bij(P_1,E);
% updating depsid
depsidm = lambdam*Aijkl_Bij(dY_dsig,df_dY);
% depsidn = lambdan*Aijkl_Bij(dY_dsig,df_dY);
% updating Omega
dOmegam = lambdam*dg_dY;
% dOmegan = lambdan*dg_dY;
Omegam = Omega+dOmegam;
% Omegan = Omega+dOmegan;
[matD,S] = matDO1(Omega);
for i=1:3
    for j=1:3
        for k=1:3
            for l=1:3
                for p=1:3
                    for q=1:3
                        dS_dOmega(i,j,k,l,p,q)=2*a1*E(i,j)*E(k,l)*E(p,q)+...
                            1/4*a2*(E(i,k)*(E(p,j)*E(q,l)+E(p,l)*E(q,j))+...
                        E(i,l)*(E(p,j)*E(q,k)+E(p,k)*E(q,j))+...
                        E(j,l)*(E(i,p)*E(q,k)+E(i,q)*E(p,k))+...
                        E(j,k)*(E(i,p)*E(q,l)+E(i,q)*E(p,l)))+...
                        1/2*a3*(E(i,j)*(E(k,p)*E(l,q)+E(k,q)*E(l,p))+E(k,l)*(E(i,p)*E(j,q)+E(i,q)*E(j,p)))+...
                        a4*(E(i,k)*E(j,l)+E(i,l)*E(j,k))*E(p,q);
                    end
                end
            end
        end
    end
end

Sm=zeros(3,3,3,3);
Sm=S+Aijklpq_Bpq(dS_dOmega,dOmegam);

sigm = sigmaT+dsig;
epsm = Aijkl_Bkl(Sm,sigm)+epsid+depsidm;
fdm = fdDP(sigm,zeros(3,3),Omegam,Pars);

end