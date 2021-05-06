% ***************************************************************************
%  Initialize all profiles and closure parameters 
%      Call once before time loop
%      Sets all forcing: pressure gradients, stresses, etc.
%      Should be used to adjust initial temperature/salinity profiles
%      Velocity initialized to zero
%      Turbulence quantities initialized to "SMALL"; Lengthscale parabolic
% ***************************************************************************

% ***************************************************************************
%  Physical parameters 
% ***************************************************************************
z0=0.01; %bottom roughness [m]
zb=10*z0; %bottom height
g=9.81; %m^2/s - gravity
C_D = 0.0025; %friction coefficient
SMALL=1e-6;
kappa=0.4; %von Karman constant
nu=1e-6; %m^2/s kinematic viscosity
u_crit = 0.007; %m/s Critical Friction velocity
Ar = 1.25; %Erosion coefficient

rhoratio = 1.01; %rhos/rhow (unitless)
D = 0.0001; %m - particle diameter

% ***************************************************************************
%  Forcing parameters - Need to modify to allow for time variable Px.
% ***************************************************************************
Px0 = -.0025; %Magnitude on pressure gradient forcing
T_Px = 0; %Periodic Pressure Forcing. Non existent in project

% ***************************************************************************
%  Turbulence closure parameters 
% ***************************************************************************
A1=0.92;
A2=0.74;
B1=16.6;
B2=10.1;
C1=0.08;
E1=1.8;
E2= 1.33;
E3=0.25;
Sq=0.2;

% ********************************************************************
%  Setup initial conditions for scalar and density
% ********************************************************************
delC=5; %change in temperature at initial themocline [deg C]; set to zero for Unstratified Case
zdelC = -5; %position of initial thermocline
dzdelC = -2; %width of initial thermocline
alpha = 0; %thermal expansivity, set to zero for passive scalar case 
rho0=1000; %kg/m^3 - water density

for i=1:N
    if i>N/2
        C(i)=0;
    else
        C(i)=0;
    end
    %if z(i)>=H/2
    %    C(i)=15;
    %elseif z(i)>=zdelC+0.5*dzdelC
    %    C(i)=15+delC;
    %else
    %    C(i)=15+delC*(z(i)-zdelC+0.5*dzdelC)/dzdelC;
    %end
   rho(i)=rho0*(1-alpha*(C(i)-15)); % Single scalar, linear equation of state
end
%Brunt-Vaisala frequency from density profile
N_BV(1)=sqrt((-g/rho0)*(rho(2)-rho(1))/(dz));
for i =2:N-1
   N_BV(i)=sqrt((-g/rho0)*(rho(i+1)-rho(i-1))/(2*dz));
end
N_BV(N)=sqrt((-g/rho0)*(rho(N)-rho(N-1))/(dz));

% *******************************************************************
%  Set initial conditions for u, q^2, and other turbulence quantities
% *******************************************************************

t(1)=0;
for i=1:N
   U(i) = 0.0;
   Up(i) = 0.0;
   V(i) = 0.0;
   Vp(i) = 0.0;
   Q2(i)=SMALL;  %"seed" the turbulent field with small values, then let it evolve
   Q2L(i)=SMALL;
   Q(i)=sqrt(Q2(i));
   L(i)=-kappa*H*(z(i)/H)*(1-(z(i)/H)); %Q2L(n,1)/Q2(n,1) = 1 at initialization;
   Gh= -(N_BV(i)*L(i)/(Q(i)+SMALL))^2; 
   Gh=min(Gh,0.0233);
   Gh=max(Gh,-0.28);
   num=B1^(-1/3)-A1*A2*Gh*((B2-3*A2)*(1-6*A1/B1)-3*C1*(B2+6*A1));
   dem=(1-3*A2*Gh*(B2+6*A1))*(1-9*A1*A2*Gh);
   Sm(i)=num/dem;
   Sh(i)=A2*(1-6*A1/B1)/(1-3*A2*Gh*(B2+6*A1)); 
   Kq(i)=Sq*Q(i)*L(i) +nu;  %turbulent diffusivity for Q2
   nu_t(i)=Sm(i)*Q(i)*L(i) +nu; %turbulent viscosity
   Kz(i)=Sh(i)*Q(i)*L(i) + nu; %turbulent scalar diffusivity
end

% *******************************************************************
%  Pre-define Tridiagonal Arrays - Just in case
% *******************************************************************
aU = zeros(1,N);
bU = zeros(1,N);
cU = zeros(1,N);
dU = zeros(1,N);
aC = zeros(1,N);
bC = zeros(1,N);
cC = zeros(1,N);
dC = zeros(1,N);
aQ2 = zeros(1,N);
bQ2 = zeros(1,N);
cQ2 = zeros(1,N);
dQ2 = zeros(1,N);
aQ2L = zeros(1,N);
bQ2L = zeros(1,N);
cQ2L = zeros(1,N);
dQ2l = zeros(1,N);


% *******************************************************************
%  Save initial conditions as first columns in saved matrix
% *******************************************************************
Um(:,1) = U;
Cm(:,1) = C;
Q2m(:,1) = Q2;
Q2Lm(:,1) = Q2L;
rhom(:,1) = rho;
Lm(:,1) = L;
nu_tm(:,1) = nu_t;
Kzm(:,1) = Kz;
Kqm(:,1) = Kq;
N_BVm(:,1) = N_BV;

% *******************************************************************
% Loading up particles
% *******************************************************************

particle_rhor = [1.05 1.5 2.25 5.0]; %m - particle density ratio
for k = 1:length(particle_rhor)
    for i = 1:pop
        part(i, k).particle = particle(particle_rhor(k));
    end
end
   
