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
h_ref = 0.01;

p = 2;
error_desired = 1e-6;
% Create butcher tableau
Bogacki = struct();
Bogacki.C = [0,1/2, 3/4, 1];
Bogacki.B = [2/9, 1/3, 4/9, 0; 7/24, 1/4, 1/3, 1/8];
Bogacki.A = [0,0,0,0; 1/2,0,0,0; 0,3/4,0,0; 2/9,1/3, 4/9, 0];

[t_list,V_list,h_avg, num_evals, num_fails] = explicit_RK_variable_step_integration(rate_func_in, t_span, V0, h_ref, Bogacki, p, error_desired);

% Calculate actual solution
solution_V_list = zeros(length(t_list), length(V0));
for i = 1:length(t_list)
    solution_V_list(i, :) = compute_planetary_motion(t_list(i), V0, orbit_params)';
end

colors = ['k', 'b', 'm', 'r'];
for i = 1:length(V0)
    plot(t_list, V_list(:,i), colors(i) + "-");
    hold on;
    plot(t_list, V_list(:,i), colors(i) + "-.");
    plot(t_list, V_list(:,i), colors(i) + ":");
    plot(t_list, V_list(:,i), colors(i) + "--");
    plot(t_list, solution_V_list(:,i), colors(i) + "-o");
end

% Step failure fraction
step_failure_rate = num_fails / (length(t_list) + num_fails - 1);

% Global truncation error (difference from exact solution)
solution_V_ad = zeros(length(t_list), length(V0));
for i = 1:length(t_list)
    solution_V_ad(i,:) = compute_planetary_motion(t_list(i), V0, orbit_params)';
end
global_error_ad = max(vecnorm(V_list - solution_V_ad, 2, 2));

fprintf('Adaptive RK:\n');
fprintf('Global error = %.3e\n', global_error_ad);
fprintf('Average step size = %.5f\n', h_avg_ad);
fprintf('Number of function evaluations = %d\n', num_evals_ad);
fprintf('Step failure rate = %.3f\n\n', step_failure_rate);

%% --- Fixed-step RK Integration ---
[t_list_fix, V_list_fix, num_evals_fix] = explicit_RK_fixed_step_integration(rate_func_in, t_span, V0, h_fixed, Bogacki);

% Global truncation error
solution_V_fix = zeros(length(t_list_fix), length(V0));
for i = 1:length(t_list_fix)
    solution_V_fix(i,:) = compute_planetary_motion(t_list_fix(i), V0, orbit_params)';
end
global_error_fix = max(vecnorm(V_list_fix - solution_V_fix, 2, 2));

fprintf('Fixed-step RK:\n');
fprintf('Global error = %.3e\n', global_error_fix);
fprintf('Average step size = %.5f\n', h_fixed);
fprintf('Number of function evaluations = %d\n', num_evals_fix);