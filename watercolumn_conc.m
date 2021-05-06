%********************************************************************
% Master code for 200B water column (vertical) cases
%    Original code from Lisa Lucas, modified by Tina Chow
%       Spring 2018: Mark Stacey and Michaella Chung
%********************************************************************

clear all;
close all;


%********************************************************************
%Define model set up - grid and timestep
%********************************************************************
pop = 10;
newpop = 25;
N=80;%number of grid points
H=20; %depth (meters)
dz=H/N; %grid spacing - may need to adjust to reduce oscillations
dt=60; %(seconds) size of time step 
M=1440; %number of time steps 
% M=1440; %number of time steps 
beta =dt/dz^2;
for i=1:N % Initialize grid   
   z(i)=-H+dz*(i-1/2); %bottom at z=-H, free surface at 0
end
isave=100; %increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
savecount=0;

% *******************************************************************
%  Call code to set up initial conditions and model parameters
% *******************************************************************
wc_setup
t = [];
% *******************************************************************
%  Start of time loop
% *******************************************************************
%tol = 0.0001;
%Uprev = U;
%wc_advance
%m = 0;
%while norm(U-Uprev,2) > tol
%Uprev = U;

figure(2);
hold on
plot(C, z);
title("Concentration Profile in Steady Flow");
xlabel('C', 'interpreter','latex');
ylabel("z");

for m=2:M
   t(m)=dt*(m-1); %de+/*ine time
   wc_advance %uses BGO/Mellor-Yamada 2-equation closure
% *******************************************************************
%  Saving Profiles - isave defines decimation
% *******************************************************************
   
   
   if mod(m,isave) == 0
      %t = [t dt*(m-1)];
      savecount = savecount+1;
      time(savecount) = isave*dt*(m-1);
      Um(:,savecount) = U;
      Cm(:,savecount) = C;
      Q2m(:,savecount) = Q2;
      Q2Lm(:,savecount) = Q2L;
      rhom(:,savecount) = rho;
      Lm(:,savecount) = L;
      nu_tm(:,savecount) = nu_t;
      Kzm(:,savecount) = Kz;
      Kqm(:,savecount) = Kq;
      N_BVm(:,savecount) = N_BV;
      
   %figure(3);
   %hold on
   %plot(C, z);
   %xlabel('C', 'interpreter','latex');
   %ylabel("z");
   %figure(4);
   %plot(U, z);
   %xlabel('C', 'interpreter','latex');
   %ylabel("z");
   end
   m = m+1;
  
end

% *******************************************************************
%  End of time loop
% *******************************************************************

% *******************************************************************
%  PLOTTING FOLLOWS HERE
%   columns of variable matrices (Um, Cm, etc) vs. z array
% *******************************************************************


figure(1);
contourf(time,z,Cm);
colorbar
title("Concentration Contours - Tidal Flow");
xlabel("t");
ylabel("z");

figure(2);
plot(Cm, z);
title("Concentration Profile in Steady Flow");
xlabel('C', 'interpreter','latex');
ylabel("z");
grid on
legend({'initial', 't1', 't2', 't3', 't4', 't5'});




% *******************************************************************
%  End of Main Program
% *******************************************************************


