% function computing the matrix Fk linking the current state with the
% following one, expressing how the state evolves in time. We assume the
% model studied during the lectures with constant velocities along x, y, z.
function [Val, Jacob] = ConstantVelocity(X, T)
    Val = zeros(size(X));
    Val(1:2:end) = X(1:2:end) + T * X(2:2:end);
    Val(2:2:end) = X(2:2:end);
    Jacob = [1,T; 0,1];
    Jacob = blkdiag(Jacob,Jacob,Jacob,Jacob);
end