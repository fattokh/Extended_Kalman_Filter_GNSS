%% Readme

% this script performs some sanity checks on kalman_filter_main to ensure that the
% code is working properly.
%
% As output we peovide:
  % - figure 1: Computed position along X, Y and Z directions
  % - figure 2: Computed velocity along X, Y and Z directions
  % - figure 3: Computed clock bias


%%

clearvars;
close all;
clc;

% number of iterations for Kalman filter and least square
iterations = 50;
dist = 2000000;

% here we set the clock bias. Start testing with this value to zero and
% then gradually increase it. The kalman filter plot for clock bias should
% converge to this value
bias = 25;

% here we set the gaussian noise added to the pseudorange measurements. Start
% testing by leaving the wgn commented and use a zero noise configuration,
% then try with the wgn (white gaussina noise). The kalman filter error
% and the Least Square error should decrease as the number of iterations
% increases

% noise = zeros(4, iterations);
noise = wgn(4, iterations, 1);

% here we add the eventual noise and clock bias to the measurements. The
% satellite positions and pseudoranges are handmade to represent the
% condition of four fixed satellites and a static receiver. This extremely
% simplified scenario is for testing purposes, to ensure that the code is
% working as intended. The receiver is supposed to be in position 0,0,0, so
% the Kalman filter and the Least Square should converge to a position
% close to 0,0,0.
for i = 1:iterations
    sat_pos(i) = {[0, 0, dist; 0, dist, 0;-dist, 0, 0; 0, -dist, 0]};
    pseudoranges(i) = {[dist+noise(1, i)+bias, dist+noise(2, i)+bias, dist+noise(3, i)+bias, dist+noise(4, i)+bias]};
end
% positioning interval
dt = 1;
f = @(X) ConstantVelocity(X, dt);

% set the covariance matrix of the gaussian noise of the state X, called Q.
% noise parameters
Sf = 36;
Sg = 0.01;
% variance of the state error w
sigma=5;

Qb = [Sf*dt+Sg*dt*dt*dt/3, Sg*dt*dt/2;
	  Sg*dt*dt/2,          Sg*dt];

Qxyz = sigma^2 * [dt^3/3, dt^2/2;
                  dt^2/2, dt];

% creates the block diagonal matrix with the variance
Q = blkdiag(Qxyz,Qxyz,Qxyz,Qb);
% initializing the state vector X    
X = zeros(8,1);
% the coordinates are estimated to a random value since this leads to a
% faster convergence of the kalman filter. Initial velocities are set to
% zero instead.
X([1 3 5]) = [50, 50, 50];
X([2 4 6]) = [0 0 0];
% initial clock bias
X(7,1) = 0;
% initial clock drift
X(8,1) = 0;
% initialization of prediction step matrix
P = eye(8)*10;
% preallocating results
pos_KF = zeros(3, iterations);
pos_LS = zeros(3, iterations);
vel_KF = zeros(3, iterations);
clock_bias = zeros(1, iterations);

% loop with the kalman filter iterations
for ii = 1:iterations
    % function to compute the pseudorange equations (looks more tidy as a function)
    g = @(X) PseudorangeEquation(X, sat_pos{ii});
    % variance of measurement error
    pseudorange_err = 36;
    % matrix of the error v
    R = eye(size(sat_pos{ii}, 1)) * pseudorange_err;
    % value of the measured pseudorange (observation)
    Z = pseudoranges{ii}.'; 
    % compute ii iteration of extended kalman filter function
    [X,P] = Extended_KF(f,g,Q,R,Z,X,P);
    % positioning using Kalman Filter
    pos_KF(:,ii) = X([1 3 5]).'; 
    % record the computed speed
    vel_KF(:,ii) = X([2 4 6]).';
    % record the computed clock bias
    clock_bias(ii) = X(7);
    % positioning using Least Square as a comparison
    pos_LS(:,ii) = Rcv_Pos_Compute(sat_pos{ii}, pseudoranges{ii});
end

% plot the estimated position x, y, z from the Kalman Filter (should converge to approx. 0,0,0)
names = ["x", "y", "z"];
for ii = 1:3
    subplot(3,1,ii)
    plot(1:iterations, pos_KF(ii,1:iterations) ,'-r')
    hold on;
    grid on;
    plot(1:iterations, pos_LS(ii,1:iterations) ,'-b')
    plot(1:iterations, zeros(iterations, 1), '-k');
    xlim([1,iterations]);
    legend('EKF','LSQ', 'True position');
    xlabel('Iterations of kalman filter')
    ylabel(names(ii)+' position [m]')
end
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf Computed position x, y and z','HorizontalAlignment','center','VerticalAlignment', 'top');

% plot the estimated velocity along the three axis, should converge to
% approximately 0,0,0 since the receiver is assumed to be fixed
figure;
for ii = 1:3
    subplot(3,1,ii)
    plot(1:iterations, vel_KF(ii,1:iterations), '-r')
    grid on;
    hold on;
    plot(1:iterations, zeros(iterations, 1), '-k');
    xlim([1,iterations]);
    legend('EKF velocity', 'True velocity');
    xlabel('Iterations of kalman filter')
    ylabel('Velocity along '+names(ii)+' [m]')
end
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf Computed velocity along x, y and z directions','HorizontalAlignment','center','VerticalAlignment', 'top');

% plot the estimated clock bias. Should be around the value set in the
% variable "bias"
figure;
plot(clock_bias);
yline(bias,'-','real value');
grid on;
xlabel('Iterations of kalman filter')
ylabel('EKF clock bias')
xlim([1,iterations]);
ylim([5,30]);

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf Computed clock bias','HorizontalAlignment','center','VerticalAlignment', 'top');
