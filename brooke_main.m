clear;
clc;
clf;
close all;

%%
orbit_params = struct();
orbit_params.m_sun = 1;
orbit_params.m_planet = 1;
orbit_params.G = 40;

rate_func_in = @(t, V) gravity_rate_func(t, V, orbit_params); 
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
[t_list, V_list, h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in, t_span, V0, h_ref, forward_euler);
[t_list, V_list_2, h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in, t_span, V0, h_ref, ralstons_method);
[t_list, V_list_3, h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in, t_span, V0, h_ref, ralstons_third_order_method);
[t_list, V_list_4, h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in, t_span, V0, h_ref, ralstons_fourth_order_method);

% Calculate actual solution
solution_V_list = zeros(length(t_list), length(V0));
for i = 1:length(t_list)
    solution_V_list(i, :) = compute_planetary_motion(t_list(i), V0, orbit_params)';
end


figure(1)
colors = ['k', 'b', 'm', 'r'];
for i = 1:length(V0)
    plot(t_list, V_list(:,i), colors(i) + "-");
    hold on;
    plot(t_list, V_list_2(:,i), colors(i) + "-.");
    plot(t_list, V_list_3(:,i), colors(i) + ":");
    plot(t_list, V_list_4(:,i), colors(i) + "--");
    plot(t_list, solution_V_list(:,i), colors(i) + "-o");
end

xlabel('Time');
ylabel('Solution');
legend("Eulers's Method", "Ralston's Method", "Ralston's 3rd-Order Method", "Ralston's 4th-Order Method", "Actual Solution")
grid on;

%% Local truncation error