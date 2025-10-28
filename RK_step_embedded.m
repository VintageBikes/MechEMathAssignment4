%This function computes the value of X at the next time step
%for any arbitrary embedded RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%OUTPUTS:
%XB1: the approximate value for X(t+h) using the first row of the Tableau
%XB2: the approximate value for X(t+h) using the second row of the Tableau
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB1, XB2, num_evals] = RK_step_embedded(rate_func_in,t,XA,h,BT_struct)
    BT_a = BT_struct.A;
    BT_b = BT_struct.B;
    BT_c = BT_struct.C;

    s = length(BT_b);
    K = zeros(length(XA), s);
    % Fill K_list
    for n = 1:s
        % sum_aK = Σ(a_i,j * K_j)
        t_temp = t + BT_c(n)*h;
        X_temp = XA + h*K*BT_a(n,:)';
        K(:, n) = rate_func_in(t_temp, X_temp);
    end
    
    % sum_bK = Σ(b_i * K_i)z
    XB1 = XA + h*K*(BT_b(1,:))';
    XB2 = XA + h*K*(BT_b(2,:))';
    num_evals = s;
end