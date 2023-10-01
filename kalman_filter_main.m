%% Readme

% Function that uses the Kalman filter for solving the GPS positioning
% problem. It makes use of the Extended_KF and least square functions.
% The input parameters are the pseudoranges and the satellites positions.
%
% We estimate the pseudorange noise v as white gaussian noise.
% The least square function can be used to calculate the unknowns x, y, z 
% and is provided as a comparison. In the KF solution we use
% the Extended Kalman filter to deal with the nonlinearity of the
% pseudorange equation. We model the system as having a constant velocity.
% The equation we are interested to solve for the GPS is
%   rho = || Xs - X || + b + v
% in which Xs and X represent the position of the satellite and
% receiver, respectively, and || Xs - X || represents the distance between 
% them. b represents the clock bias of receiver, and it need to be solved 
% along with the position of receiver. rho is a measurement given by 
% receiver for each satellites, and v is the pseudorange measurement noise 
% modeled as white noise.
%
% As output we peovide:
%   - figure 1: Comparison between the error of the Kalman Filter and the
%               Least square when used to compute the position.
%
% References:
% 1. R G Brown, P Y C Hwang, "Introduction to random signals and applied 
%   Kalman filtering : with MATLAB exercises and solutions",1996
% 2. Pratap Misra, Per Enge, "Global Positioning System Signals, 
%   Measurements, and Performance(Second Edition)",2006

%% Reading the initial data

clearvars;
close all;
clc;

% position of satellites: each cell is a 2D matrix containing the x,y,z position of each
% one of the 4 satellites.
load sat_pos; 
% pseudoranges: each cell has the 4 pseudorange measurements from the 4 satellites
load pseudoranges 

% positioning interval (interval between new messages by the communication team)
dt = 1;
% kalman filter iterations (number of data slot received by previous group)
iterations = 25;

% the vector of the state to estimate is as [x vx y vy z vz b d] column
% vector, where x, y, z are the coordinates of the position, and vx, vy, vz
% the velocities on the three dimensions, b and d represent the clock bias
% and its drift rate.

f = @(X) ConstantVelocity(X, dt);

% set the covariance matrix of the gaussian noise of the state X, called Q.
% noise parameters
Sf = 36;
Sg = 0.01;
% variance of the state error w
sigma=1;

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
X([1 3 5]) = [-2.168816181271560e+006, 4.386648549091666e+006, 4.077161596428751e+006];
X([2 4 6]) = [0 0 0];
% initial clock bias [microseconds]
X(7,1) = 3.575261153706439e+006;
% initial clock drift
X(8,1) = 4.549246345845814e+001;
% initialization of prediction step matrix
P = eye(8)*10;

pos_KF = zeros(3, iterations);
pos_LS = zeros(3, iterations);
% loop with the kalman filter iterations
for ii = 1:iterations
    % function to compute the pseudorange equations (looks more tidy as a function)
    g = @(X) PseudorangeEquation(X, SV_Pos{ii});
    % variance of measurement error
    pseudorange_err = 36;
    % matrix of the error v
    R = eye(size(SV_Pos{ii}, 1)) * pseudorange_err;
    % value of the measured pseudorange (observation)
    Z = SV_Rho{ii}.'; 
    % compute ii iteration of extended kalman filter function
    [X,P] = Extended_KF(f,g,Q,R,Z,X,P);
    % positioning using Kalman Filter
    pos_KF(:,ii) = X([1 3 5]).'; 
    % positioning using Least Square as a contrast
    pos_LS(:,ii) = Rcv_Pos_Compute(SV_Pos{ii}, SV_Rho{ii});
end

%% Plotting the results

% plot of the position error for both KF and LS
% in order to plot the relative error it is sufficient to substract the mean
% from the position estimation, since we do not precisely know the moon
% surface ground.
names = ["x", "y", "z"];
for ii = 1:3
    subplot(3,1,ii)
    plot(1:iterations, abs(pos_KF(ii,1:iterations) - mean(pos_KF(ii,1:iterations))),'-r')
    hold on;
    grid on;
    plot(1:iterations, abs(pos_LS(ii,1:iterations) - mean(pos_KF(ii,1:iterations))),'-b')
    legend('EKF','LSQ')
    xlabel('Iterations')
    ylabel('Error along '+names(ii)+' [m]')
end
ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,'\bf Relative positioning error in x,y and z directions','HorizontalAlignment','center','VerticalAlignment', 'top');
