% function to compute least square position for GPS.
% @input sat_pos:  position of satellites
% @input pseudoranges: pseudorange of satellites
% @output receiver_pos: position of receiver 
% @output clk_bias: receiver clock bias

function [receiver_pos,clk_bias] = Rcv_Pos_Compute(sat_pos, pseudoranges)
    % number of visible satellites
    num_sats = size(sat_pos);
    % initialize GPS position guess and clock bias
    receiver_pos=[0 0 0];
    clk_bias=0;
    % checking if there are always at least four available satellites 
    if num_sats<4
        return;
    end
    % start iterating
    % constraint for convergence  
    flag=1;
    exit_loop=100;
    count=0;
    
    % the iterative approch will go on until convergence is achieved, i.e. the deviation is 
    % the smallest possible. Which is limited by the accuracy of the system of measurement.   
    while (exit_loop > flag) 
        % compute G, the Jacobian unitary matrix
        G = G_Compute(sat_pos, receiver_pos);
        % compute delta_rho (linearization of pseudorange)
        delta_rho = Delta_Rho_Compute(pseudoranges, sat_pos, receiver_pos, clk_bias);
        % iteration for new postion (solution)
        delta_X = inv(G' * G) * G' * delta_rho;
        receiver_pos = (receiver_pos' + delta_X(1:3))';
        clk_bias = clk_bias + delta_X(4);
        % update exit condition
        exit_loop = (delta_X(1)^2 + delta_X(2)^2 + delta_X(3)^2)^0.5;
        count = count+1;
        % if after ten iterations we do not exit from the loop, exit anyway
        if count>10
            exit_loop=flag/2;
        end
    end
end


% function to compute the Jacobian unitary matrix.
% @input SV_Pos: position of satellites
% @input Rcv_Pos: position of receiver 
% @output : G = [ a_x(1)   a_y(1)   a_z(1)   1 
%                 a_x(2)   a_y(2)   a_z(2)   1 
%                 a_x(3)   a_y(3)   a_z(3)   1 
%                 a_x(4)   a_y(4)   a_z()4   1 ]
%
% the vector a = [ a_x(i)   a_y(i)   a_z(i) ] is the unitary vector connecting the
% user to the i-th satellite

function G = G_Compute(SV_Pos, Rcv_Pos)
    % m = number of visible satellites 
    [m, ~] = size(SV_Pos);
    % element-wise difference between the arrays
    dX = bsxfun(@minus, SV_Pos, Rcv_Pos);
    % norm
    nor = sum(dX .^2, 2) .^0.5;
    % normalization 
    unit_matrix = bsxfun(@rdivide, dX, nor);
    G = [-unit_matrix ones(m,1)];
end


% function to compute the delta_rho
% delta_rho comes from the linearization of pseudorange equations 
% @input pseudoranges: pseudorange of satellites
% @input sat_pos:  position of satellites
% @input rec_pos: receiver position
% @input clk_bias: receiver clock bias 
% @output : delta_rho 

function [delta_rho] = Delta_Rho_Compute(pseudoranges, sat_pos, rec_pos, clk_bias)
    % m = number of visible satellites
    [m, ~] = size(sat_pos);
    delta_rho = zeros(m, 1);
    for i = 1:m
         % pseudoragne estimation 
        rho0 = norm(sat_pos(i,:) - rec_pos) + clk_bias;
        delta_rho(i,1) = pseudoranges(i) - rho0;
    end
end