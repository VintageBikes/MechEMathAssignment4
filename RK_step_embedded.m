function [XB1, XB2, num_evals] = RK_step_embedded(rate_func_in, t, XA, h, BT_struct)
    BT_a = BT_struct.A;
    BT_b = BT_struct.B;
    BT_c = BT_struct.C;

    s = length(BT_c);
    K = zeros(length(XA), s);
    num_evals = 0;

    % Fill K_list
    for n = 1:s
        t_temp = t + BT_c(n)*h;
        X_temp = XA + h*K*BT_a(n,:)';
        K(:, n) = rate_func_in(t_temp, X_temp);
        num_evals = num_evals + 1;
    end
    
    % sum_bK = Σ(b_i * K_i)
    XB1 = XA + h*K*(BT_b(1,:))';
    XB2 = XA + h*K*(BT_b(2,:))';
end
