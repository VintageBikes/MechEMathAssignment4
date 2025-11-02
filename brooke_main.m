clear;
clc;
clf;
close all;

%% Calculate RK Solutions and Actual Solution
orbit_params = struct();
orbit_params.m_sun = 1;
orbit_params.m_planet = 1;
orbit_params.G = 40;

rate_func_in = @(t, V) gravity_rate_func(t, V, orbit_params);
solution_func_in = @(t, V) compute_planetary_motion(t, V, orbit_params);
t_span = [0, 16];
V0 = [8; 0; 0; 1.5];
h_ref = 0.2;

% Create butcher tableau
% euler's
forward_euler = struct();
forward_euler.A = [0];
forward_euler.B = [1];
forward_euler.C = [0];
% ralston's
ralstons_method = struct();
ralstons_method.A = [0 0; 2/3 0];
ralstons_method.B = [1/4 3/4];
ralstons_method.C = [0; 2/3];
% ralston's third order
ralstons_third_order_method = struct();
ralstons_third_order_method.A = [0 0 0; 1/2 0 0; 0 3/4 0];
ralstons_third_order_method.B = [2/9 1/3 4/9];
ralstons_third_order_method.C = [0; 1/2; 3/4];
% ralston's fourth order
ralstons_fourth_order_method = struct();
ralstons_fourth_order_method.A = [
    0 0 0 0
    2/5 0 0 0
    (-2889+1428*sqrt(5))/1024 (3785-1620*sqrt(5))/1024 0 0
    (-3365+2094*sqrt(5))/6040 (-975-3046*sqrt(5))/2552 (467040+203968*sqrt(5))/240845 0
];
ralstons_fourth_order_method.B = [
    (263+24*sqrt(5))/1812 (125-1000*sqrt(5))/3828 (3426304+1661952*sqrt(5))/5924787 (30-4*sqrt(5))/123
];
ralstons_fourth_order_method.C = [0 2/5 (14-3*sqrt(5))/16 1];

% Calculate RK solutions
[t_list, V_list_1, h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in, t_span, V0, h_ref, forward_euler);
[t_list, V_list_2, h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in, t_span, V0, h_ref, ralstons_method);
[t_list, V_list_3, h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in, t_span, V0, h_ref, ralstons_third_order_method);
[t_list, V_list_4, h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in, t_span, V0, h_ref, ralstons_fourth_order_method);

% Calculate actual solution
solution_V_list = zeros(length(t_list), length(V0));
for i = 1:length(t_list)
    solution_V_list(i, :) = compute_planetary_motion(t_list(i), V0, orbit_params)';
end

%% Approximate Solutions vs. Actual Solutions
figure(1)
colors = ['k', 'b', 'm', 'r'];
for i = 1:length(V0)
    plot(t_list, V_list_1(:,i), colors(i) + "-");
    hold on;
    plot(t_list, V_list_2(:,i), colors(i) + "-.");
    plot(t_list, V_list_3(:,i), colors(i) + ":");
    plot(t_list, V_list_4(:,i), colors(i) + "--");
    plot(t_list, solution_V_list(:,i), colors(i) + "-o");
end

title("Runge-Kutta Methods vs. Time")
xlabel('Time');
ylabel('Solution');
legend("Eulers's Method", "Ralston's Method", "Ralston's 3rd-Order Method", "Ralston's 4th-Order Method", "Actual Solution")
grid on;

%% Calculate Local Truncation Error
num_points = 100;
t_ref = 0.492;
h_list = logspace(-5, 1, num_points);

te_1 = local_truncation_error(forward_euler, rate_func_in, solution_func_in, V0, t_ref, h_list);
te_2 = local_truncation_error(ralstons_method, rate_func_in, solution_func_in, V0, t_ref, h_list);
te_3 = local_truncation_error(ralstons_third_order_method, rate_func_in, solution_func_in, V0, t_ref, h_list);
te_4 = local_truncation_error(ralstons_fourth_order_method, rate_func_in, solution_func_in, V0, t_ref, h_list);

% Calculate analytical difference
analytical_difference = zeros(length(h_list), 1);
for i = 1:length(h_list)
    analytical_difference(i) = norm(solution_func_in(t_ref + h_list(i), V0) - solution_func_in(t_ref, V0));
end

% Create Fit Lines
p_analytical = polyfit(log(h_list), log(analytical_difference), 1);
y_hat_analytical = exp(p_analytical(1) * log(h_list) + p_analytical(2));

threshold = 4e-15;  % h_values are stable until they exceed this threshold
ix = find(h_list >= 5, 1);  % ignore h_values greater than 5 for fit lines

iy_1 = find(te_1 > threshold, 1);
p_1 = polyfit(log(h_list(iy_1:ix)), log(te_1(iy_1:ix)), 1);
y_hat_1 = exp(p_1(1) * log(h_list) + p_1(2));

iy_2 = find(te_2 > threshold, 1);
p_2 = polyfit(log(h_list(iy_2:ix)), log(te_2(iy_2:ix)), 1);
y_hat_2 = exp(p_2(1) * log(h_list) + p_2(2));

iy_3 = find(te_3 > threshold, 1);
p_3 = polyfit(log(h_list(iy_3:ix)), log(te_3(iy_3:ix)), 1);
y_hat_3 = exp(p_3(1) * log(h_list) + p_3(2));

iy_4 = find(te_4 > threshold, 1);
p_4 = polyfit(log(h_list(iy_4:ix)), log(te_4(iy_4:ix)), 1);
y_hat_4 = exp(p_4(1) * log(h_list) + p_4(2));

% Plot
figure(2);
loglog(h_list, y_hat_analytical,'-k')
hold on;
loglog(h_list, y_hat_1, '-b')
loglog(h_list, y_hat_2, '-m')
loglog(h_list, y_hat_3, '-r')
loglog(h_list, y_hat_4, '-g')

loglog(h_list, analytical_difference, '+k')
loglog(h_list, te_1, '+b');
loglog(h_list, te_2, '+m');
loglog(h_list, te_3, '+r');
loglog(h_list, te_4, '+g');
xlabel('Step Size (h)'); ylabel('Local Truncation Error'); title('Local Truncation Error for Runge-Kutta'); 
legend("Analytical Difference", "Eulers's Method", "Ralston's Method", "Ralston's 3rd-Order Method", "Ralston's 4th-Order Method");
grid on;

%% Local Truncation Error Table



%% Local Truncation Error Function
function truncation_error = local_truncation_error(BT_struct, test_func, solution_func, V0, t_ref, h_list)    
    truncation_error = zeros(length(h_list), 1);
    for i = 1:length(h_list)
        X_next = explicit_RK_step(test_func, t_ref, solution_func(t_ref, V0), h_list(i), BT_struct);
        X_guess = solution_func(t_ref + h_list(i), V0);
        truncation_error(i) = norm(X_next-X_guess);
    end
end