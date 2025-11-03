function day_2025_10_21_coding()
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;

    X0 = 1;

    my_rate = @(t_in, X_in) rate_func01(t_in,X_in);
    tspan = [0,30];

    DormandPrince = struct();
    DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
    5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DormandPrince.A = [0,0,0,0,0,0,0;
    1/5, 0, 0, 0,0,0,0;...
    3/40, 9/40, 0, 0, 0, 0,0;...
    44/45, -56/15, 32/9, 0, 0, 0,0;...
    19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
    9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
    35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];

    
    h_ref = 0.5;
    
    %[t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(my_rate,tspan,V0,h_ref,DormandPrince);
    
    
    % figure(1)
    % subplot(2, 1, 1);
    % hold on;
    % plot(t_range,V_list(:, 1), 'k', 'linewidth', 2);
    % plot(t_range,V_list(:, 2), 'b', 'linewidth', 2);
    % 
    % plot(t_list, X_list(:, 1), 'r--', 'linewidth', 2);
    % plot(t_list, X_list(:, 2), 'r--', 'linewidth', 2);
    % xlabel('time');
    % ylabel('position component');
    % 
    % subplot(2, 1, 2);
    % hold on;
    % plot(t_range,V_list(:, 3), 'k', 'linewidth', 2);
    % plot(t_range,V_list(:, 4), 'b', 'linewidth', 2);
    % 
    % plot(t_list, X_list(:, 3), 'r--', 'linewidth', 2);
    % plot(t_list, X_list(:, 4), 'r--', 'linewidth', 2);
    % xlabel('time');
    % ylabel('velocity component');


    %local truncation error experiments for embedded method
    n_samples = 60;
    h_ref_list = logspace(-3, 1, n_samples);

    error_proxy_list = zeros(1, n_samples);

     for n = 1:length(h_ref_list)
         h_ref = h_ref_list(n);
         
         [XB1,XB2,~] = RK_step_embedded(my_rate, tspan(1), X0, h_ref, DormandPrince);
            
         error_proxy_list(n) = abs(XB1 - XB2);
     end

    filter_params = struct();
    filter_params.min_yval = 1e-13;
    filter_params.max_yval = 1e-6;

    [p_proxy, k_proxy] = loglog_fit(h_ref_list, error_proxy_list, filter_params);


    figure(2);
    loglog(h_ref_list,error_proxy_list,'mo','MarkerFaceColor', 'm', 'MarkerSize',2);
    hold on
    loglog(h_ref_list, k_proxy * h_ref_list.^p_proxy, 'm','LineWidth', 1);
    xlabel("Step Size h");
    ylabel("|XB1 - XB2|");
    title('Embedded RK Error Proxy vs Step Size');
    legend('Error Proxy')
hold off

end

function dXdt = rate_func01(t,X)
    dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
    X = cos(t);
end