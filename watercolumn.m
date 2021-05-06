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
pop = 100;
newpop = 0;
N=80;%number of grid points
H=5; %depth (meters)
dz=H/N; %grid spacing - may need to adjust to reduce oscillations
dt=60; %(seconds) size of time step 
M=900; %number of time steps 
% M=1440; %number of time steps 
beta =dt/dz^2;
for i=1:N % Initialize grid   
   z(i)=-H+dz*(i-1/2); %bottom at z=-H, free surface at 0
end
isave=10; %increments for saving profiles. set to 1 to save all; 10 saves every 10th, etc. 
savecount=0;

% *******************************************************************
%  Call code to set up initial conditions and model parameters
% *******************************************************************
wc_setup

disp('Please enter which case you are working with:');
inp = input('1) Concentration-Based Model \n2) Particle Based \nInput: ');

% *******************************************************************
%  Start of time loop
% *******************************************************************

filename = 'projectTest.gif';
nImages = 0;
meanU = [];

if inp==2
    fignum=3;
    for m=2:M
       t(m)=dt*(m-1); %define time
       wc_advance %uses BGO/Mellor-Yamada 2-equation closure
    % *******************************************************************
    %  Saving Profiles - isave defines decimation
    % *******************************************************************

    % *******************************************************************
    %  Plotting particles
    % *******************************************************************   
        for k = 1:size(part,2)
            for i = 1:size(part,1)
               partx(i, k) = part(i, k).particle.x;
               partz(i, k) = part(i, k).particle.z;
            end
        end
        
       %particle_plot_gif(part, partx, partz, U, meanU, H, m, dt, z, nImages)
        
       if mod(m,isave) == 0
           figure(fignum)
            particle_snapshot(part, partx, partz, U, H, m, dt)
           fignum=fignum+1;
       end
    end
    
    particle_snapshot(part, partx, partz, U, H, m, dt)
    fignum=fignum+1;

    % Saving animation
    %for idx = 1:nImages
    %    [A, map] = rgb2ind(im{idx}, 256);
    %    if idx == 1
    %        imwrite(A, map, filename, 'gif', 'LoopCount', Inf,'DelayTime', 1);
    %    else 
    %        imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    %    end
    %end
    
    %figure(2);
    %t_hr = t./3600;
    %plot(t_hr(2:end), meanU);
    %title("Flow Velocity vs Time");
    %xlabel('t (hours)', 'interpreter','latex');
    %ylabel("U (m/s)");
    %grid on

% *******************************************************************
%  End of Main Program
% *******************************************************************
elseif inp==1
    figure(2);
    hold on
    plot(C, z);
    title("Concentration Profile");
    xlabel('C', 'interpreter','latex');
    ylabel("z");
    grid on
    for m=2:M
       t(m)=dt*(m-1); %de+/*ine time
       wc_advance %uses BGO/Mellor-Yamada 2-equation closure
    % *******************************************************************
    %  Saving Profiles - isave defines decimation
    % *******************************************************************
        figure(4)
        subplot(2,1,1);
        plot(C, z);
        title("Concentration Profile");
        xlabel('C', 'interpreter','latex');
        ylabel("z");

        subplot(2,1,2);
        plot(U,z);
        ylim([-H 0])
        xlim([min(U) max(U)])
        xlabel('U');
        ylabel('z');
        title('Flow Velocity');

       if mod(m,isave) == 0
          %t = [t dt*(m-1)];
          savecount = savecount+1;
          time(savecount) = (dt*(m))/3600;
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
       end
       m = m+1;
    end
    figure(1);
    contourf(time,z,Cm);
    colorbar
    title("Concentration Profile in Shallow Channel - Low Density");
    xlabel("t (hour)");
    ylabel("z (m)");

    figure(2);
    plot(Cm, z);
    title("Concentration Profile in Shallow Channel - Low Density");
    xlabel('C', 'interpreter','latex');
    ylabel("z");
    grid on
    legend({'initial', 't1', 't2', 't3', 't4', 't5'});
end

% *******************************************************************
%  End of time loop
% *******************************************************************

% *******************************************************************
%  PLOTTING FOLLOWS HERE
%   columns of variable matrices (Um, Cm, etc) vs. z array
% *******************************************************************

% Particle plot function

function particle_plot_gif(part, partx, partz, U, meanU, H, m, dt, z, nImages)
        figure(1)
        subplot(2,1,1);
        c = ['r';'k';'g';'b'];
        sh = ['+'; 'o'; 'x'; '*']; 
        for k = 1:size(part,2)
            s = scatter(partx(:,k), partz(:,k), 3, c(k));
            s.Marker = sh(k);
            hold on
        end
        hold off
        l = legend({'$\rho_{s} / \rho_{w}$=1.05', '$\rho_{s} / \rho_{w}$=1.5', '$\rho_{s} / \rho_{w}$=2.25', '$\rho_{s} / \rho_{w}$=5.0'});
        set(l, 'Interpreter', 'Latex');
        ylim([-H 0])
        xlim([0 60000])
        steps = strcat(num2str((m*dt)/3600), ' h');
        xlabel('x');
        ylabel('z');
        text(15000,-H/2, steps);
        grid on
        title('Particle Tracking - Reservoir Entrance Flow');
        subplot(2,1,2);
        plot(U,z);
        ylim([-H 0])
        xlim([min(U) max(U)])
        xlabel('U');
        ylabel('z');
        title('Flow Velocity');
        meanU = [meanU mean(U)];

        nImages = nImages+1;
        drawnow
        frame = getframe(figure(1));
        im{nImages} = frame2im(frame);
end

function particle_snapshot(part, partx, partz, U, H, m, dt)
        
         list = ['m';'k';'r';'b'];
         sh = ['+'; 'o'; 'x'; '*']; 
        for k = 1:size(part,2)
            v = scatter(partx(:,k), partz(:,k), 25, list(k));
            v.Marker = sh(k);
            hold on
        end
        hold off
        l = legend({'$\rho_{s} / \rho_{w}$=1.05', '$\rho_{s} / \rho_{w}$=1.5', '$\rho_{s} / \rho_{w}$=2.25', '$\rho_{s} / \rho_{w}$=5.0'});
        set(l, 'Interpreter', 'Latex');
        ylim([-H 0])
        xlim([0 60000])
        steps = strcat(num2str((m*dt)/3600), ' h');
        veloce = strcat(num2str(mean(U)), ' m/s');
        xlabel('x');
        ylabel('z');
        text(15000,-H/2, steps);
        text(15000,(-H/2)-0.5, veloce);
        grid on
        title('Particle Tracking - Reservoir Entrance Flow');
       
end
