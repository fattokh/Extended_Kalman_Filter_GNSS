% function to compute the GPS equation rho = || Xs - X || + b and its Jacobian.
function [Val, Jacob] = PseudorangeEquation(X, SV)
    % this executes the subtraction between X (user position) and Xs taken from
    % SV containing the satellite position. It results in the user-satellite
    % distance
    dX = bsxfun(@minus, X([1,3,5])', SV);% X - Xs
    % norm in three dimensions: sqrt((x-xsat)^2+(y-ysat)^2+...) plus time drift
    Val = sum(dX .^2, 2) .^0.5 + X(7);
    % compute the jacobian matrix for the filter
    Jacob = zeros(size(SV, 1), size(X, 1));
    % divides all the elements of dX for val
    Jacob(:, [1,3,5]) = bsxfun(@rdivide, dX, Val);
    % constant time drift
    Jacob(:, 7) = 1;
end