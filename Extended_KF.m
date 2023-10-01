% function implementing the extended Kalman filter
%
% state equation: X(n+1) = f(X(n)) + w(n+1) as seen in the lecture
%
% observation equation:
%     Z(n) = h(X(n)) + v(n)
%   where:
%     w ~ N(0,Q) gaussian noise with average null and covariance matrix Q
%     v ~ N(0,R) gaussian noise with average null and covariance matrix R
%
% input:
%     f: function for state transition, it takes a state variable Xn and
%       returns f(Xn) and the Jacobian of f at Xn
%     h: function for measurement, it takes the state variable Xn and
%       returns g(Xn) and the Jacobian of g at Xn.
%     Q: process noise covariance matrix
%     R: measurement noise covariance matrix
%     Z: current measurement
%     Xi: "a priori" state estimate (the current estimation of the state)
%     Pi: "a priori" estimated state covariance
%
% output:
%     Xo: "a posteriori" state estimate (the next estimation of the state)
%     Po: "a posteriori" estimated state covariance
%
% algorithm for extended Kalman filter:
% linearize input functions f and g to get fy (state transition matrix)
% and H(observation matrix), then apply an ordinary Kalman Filter:
%
% state equation:
%     X(n+1) = fy * X(n) + w(n), fy is the output of the linearized
%       update function.
% 
% observation equation:
%     Z(n) = H * X(n) + v(n)
%
% 1. Xp = f(Xi)               : one step projection, provides 
%                               linearization point Xp
% 2. 
%       d f |
% fy = -----|                 : linearize state equation, fy is the
%       d X |X=Xp               Jacobian of the process model
% 3.
%       d h |
% H  = -----|                 : linearize observation equation, H is
%       d X |X=Xp               the Jacobian of the measurement model
% 4. Pp = fy * Pi * fy' + Q   : covariance of Xp
% 5. K = Pp * H' * inv(H * Pp * H' + R): Kalman gain
% 6. Xo = Xp + K * (Z - g(Xp)): output state
% 7. Po = [I - K * H] * Pp    : covariance of Xo
	                                                                            
function [Xo,Po] = Extended_KF(f,h,Q,R,Z,Xi,Pi)
    N_state = size(Xi, 1);    
    [Xp, ~] = f(Xi); % step 1
    [~, fy] = f(Xp); % step 2
    [hXp, H] = h(Xp); % step 3
    Pp = fy * Pi * fy.' + Q; % step 4
    K = Pp * H' / (H * Pp * H.' + R); % step 5
    Xo = Xp + K * (Z - hXp); % step 6
    I = eye(N_state, N_state);
    Po = (I - K * H) * Pp; % step 7
end
 