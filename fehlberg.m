function day_2025_10_21_coding()
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;

    V0 = [x0;y0;dxdt0;dydt0];
    tspan = [0,30];
    t_range = linspace(tspan(1), tspan(2), 100);
    V_list = compute_planetary_motion(t_range, V0, orbit_params);

    my_rate = @(t_in, V_in) gravity_rate_func(t_in, V_in, orbit_params);

    Fehlberg = struct();
    Fehlberg.C = [0, 1/4, 3/8, 12/13, 1, 1/2];
    Fehlberg.B = [16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55;...
    25/216, 0, 1408/2565, 2197/4104, -1/5, 0];
    Fehlberg.A = [0,0,0,0,0,0;...
    1/4, 0,0,0,0,0;...
    3/32, 9/32, 0,0,0,0;...
    1932/2197, -7200/2197, 7296/2197, 0,0,0;...
    439/216, -8, 3680/513, -845/4104, 0,0;...
    -8/27, 2, -3544/2565, 1859/4104, -11/40, 0];

     h_ref = 0.5;
    
  

    %local truncation error experiments for embedded method
    n_samples = 60;
    h_ref_list = logspace(-3, 1, n_samples);

    abs_diff_list = zeros(1, n_samples);
    tr_error_list1 = zeros(1, n_samples);
    tr_error_list2 = zeros(1, n_samples);

     for n = 1:length(h_ref_list)
         h_ref = h_ref_list(n);
         V_list = compute_planetary_motion(tspan(1)+h_ref, V0, orbit_params);

         [XB1,XB2,~] = RK_step_embedded(my_rate, tspan(1), V0, h_ref, Fehlberg);

         abs_diff_list(n) = norm(V_list-V0);
         tr_error_list1(n) = norm(XB1-V_list);
         tr_error_list2(n) = norm(XB2-V_list);
     end

    filter_params = struct();
    filter_params.min_yval = 1e-13;
    filter_params.max_yval = 1e-6;

    [p1, k1] = loglog_fit(h_ref_list, tr_error_list1, filter_params);
    [p2, k2] = loglog_fit(h_ref_list, tr_error_list2, filter_params);

    p1;
    p2;

    figure(2);
    loglog(h_ref_list,abs_diff_list , 'ro', 'MarkerFaceColor', 'r','MarkerSize',2);
    hold on
    loglog(h_ref_list,tr_error_list1 , 'bo', 'MarkerFaceColor', 'b','MarkerSize',2);
    loglog(h_ref_list,tr_error_list2 , 'go', 'MarkerFaceColor', 'g','MarkerSize',2);

    loglog(h_ref_list, k1*h_ref_list.^p1, 'k', 'linewidth', 1);
    loglog(h_ref_list, k2*h_ref_list.^p2, 'k', 'linewidth', 1);
hold off

end
