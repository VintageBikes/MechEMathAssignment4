clear;
clc;
clf;
close all;

%% Calculate RK Solutions and Actual Solution
orbit_params = struct();
orbit_params.m_sun = 1;
orbit_params.m_planet = 1/330000;
orbit_params.G = 40;

rate_func_in = @(t, V) gravity_rate_func(t, V, orbit_params);
solution_func_in = @(t, V) compute_planetary_motion(t, V, orbit_params);
t_span = [0, 16];
V0 = [0.5; 0; -5; 10.28];
h_ref = 0.03;

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


figure()
plot(V_list(:, 1), V_list(:, 2), 'b-')
hold on
plot(solution_V_list(:, 1), solution_V_list(:, 2), 'k-')
plot(0,0, 'ro', 'LineWidth',5)
legend('Variable Step Approx', 'Actual Solution')
xlim([-2.5, 2.5])
ylim([-2.5, 2.5])
title('Variable Step Elliptical Orbit')
grid on;

%% --- Fixed-step RK (Ralston 3rd-order) ---
BT_fixed = struct();
BT_fixed.A = [0 0 0;
              1/2 0 0;
              0 3/4 0];
BT_fixed.B = [2/9, 1/3, 4/9];
BT_fixed.C = [0; 1/2; 3/4];

[t_fix, V_fix, h_avg_fix, num_evals_fix] = explicit_RK_fixed_step_integration(rate_func_in, t_span, V0, h_ref, BT_fixed);

%% --- Variable-step RK (Bogacki-Shampine) ---
BT_var = struct();
BT_var.A = [0 0 0 0;
              1/2 0 0 0;
              0 3/4 0 0;
              2/9 1/3 4/9 0];
BT_var.B = [2/9, 1/3, 4/9, 0;      % main solution
              7/24, 1/4, 1/3, 1/8];   % embedded solution for error
BT_var.C = [0; 1/2; 3/4; 1];

[t_var, V_var, h_avg_var, num_evals_var, num_fails] = ...
    explicit_RK_variable_step_integration(rate_func_in, t_span, V0, h_ref, BT_var, p, error_desired);

%% --- Plotting Orbits ---
% Actual solution
solution_V_list = zeros(length(t_var), length(V0));
for i = 1:length(t_var)
    solution_V_list(i,:) = compute_planetary_motion(t_var(i), V0, orbit_params)';
end

figure()
plot(V_var(:,1), V_var(:,2), 'b-', 'LineWidth', 1.5); hold on;
plot(V_fix(:,1), V_fix(:,2), 'r--', 'LineWidth', 1.5);
plot(solution_V_list(:,1), solution_V_list(:,2), 'k-', 'LineWidth', 1.5);
plot(0,0,'ro','LineWidth',3); % sun
legend('Variable-step RK','Fixed-step RK','Analytical','Sun');
xlim([-2.5,2.5]); ylim([-2.5,2.5]);
title('Planetary Motion Orbits');
xlabel('x'); ylabel('y'); grid on; axis equal;

% %% --- Display stats ---
% fprintf('Fixed-step RK:\n');
% fprintf('Average step size = %.5f\n', h_avg_fix);
% fprintf('Number of function evaluations = %d\n\n', num_evals_fix);
% 
% fprintf('Variable-step RK:\n');
% fprintf('Average step size = %.5f\n', h_avg_adapt);
% fprintf('Number of function evaluations = %d\n', num_evals_adapt);
% fprintf('Failed Step Rate = %d\n', num_fails);
%% --- Global Truncation Error Analysis (Fixed vs Adaptive) ---
h_list = logspace(-3, -1, 100);  % reference step sizes
gte_fixed = zeros(length(h_list),1);
gte_var = zeros(length(h_list),1);
h_avg_fixed_list = zeros(length(h_list),1);
h_avg_var_list = zeros(length(h_list),1);
step_failure_rate_list = zeros(length(h_list),1);

for i = 1:length(h_list)
    %Fixed step
    [t_fix, V_fix, h_avg_fix, num_evals_fix] = explicit_RK_fixed_step_integration(rate_func_in, t_span, V0, h_list(i), BT_fixed);
    v_actual = compute_planetary_motion(t_fix(end), V0, orbit_params);
    gte_fixed(i) = norm(V_fix(end,:)' - v_actual);
    h_avg_fixed_list(i) = h_avg_fix;
    
    %Variable step
    [t_var, V_var, h_avg_var, num_evals_var, num_fails] = ...
        explicit_RK_variable_step_integration(rate_func_in, t_span, V0, h_list(i), BT_var, p, error_desired);
    v_actual = compute_planetary_motion(t_var(end), V0, orbit_params);
    gte_var(i) = norm(V_var(end,:)' - v_actual);
    h_avg_var_list(i) = h_avg_var;
    
    % Step failure fraction
    step_failure_rate_list(i) = num_fails;
end

%% --- Plot Global Truncation Error ---
figure();
loglog(h_avg_fixed_list, gte_fixed, 'b*'); hold on;
loglog(h_avg_var_list, gte_var, 'r*');
xlabel('Average Step Size (h)');
ylabel('Global Truncation Error');
title('Global Truncation Error: Fixed vs Variable Step');
legend('Fixed-step RK','Variable-step RK','Location','best');
grid on;

%% --- Global Truncation Error vs Number of Function Evaluations ---
num_evals_fixed_list = zeros(length(h_list),1);
num_evals_var_list   = zeros(length(h_list),1);

for i = 1:length(h_list)
    % Fixed-step
    [~, ~, ~, num_evals_fix] = explicit_RK_fixed_step_integration(rate_func_in, t_span, V0, h_list(i), BT_fixed);
    num_evals_fixed_list(i) = num_evals_fix;

    % Variable-step
    [~, ~, ~, num_evals_var, ~] = explicit_RK_variable_step_integration(rate_func_in, t_span, V0, h_list(i), BT_var, p, error_desired);
    num_evals_var_list(i) = num_evals_var;
end

% Plot
figure();
loglog(num_evals_fixed_list, gte_fixed, 'b*'); hold on;
loglog(num_evals_var_list, gte_var, 'r*');
xlabel('Number of Function Evaluations');
ylabel('Global Truncation Error');
title('Global Truncation Error vs Function Evaluations');
legend('Fixed Step','Variable Step');
grid on;


%YURRRR
%%AVG Step vs Step Fail Rate
figure();
semilogx(h_avg_var_list, step_failure_rate_list, 'r*');
xlabel('Average Step Size (h)');
ylabel('Step Failure Count');
title('Step Failure vs Average Step Size (Variable-step RK)');
grid on;

%%Position vs time and Velocity vs tuime
figure
subplot(2,1,1);
plot(t_list, V_list(:,1),'b-','markerfacecolor','k','markeredgecolor','k','markersize',2);
hold on;
plot(t_list, V_list(:,2),'r-','markerfacecolor','k','markeredgecolor','k','markersize',2);
plot(t_list, V_list(:,1),'k.','LineWidth', 1);
plot(t_list, V_list(:,2),'k.','LineWidth', 1);
legend('X', 'Y')
title('Position')

subplot(2,1,2);
plot(t_list, V_list(:,3),'b-','markerfacecolor','k','markeredgecolor','k','markersize',10);
hold on;
plot(t_list, V_list(:,4),'r-','markerfacecolor','k','markeredgecolor','k','markersize',2);
plot(t_list, V_list(:,3),'k.','LineWidth', 1);
plot(t_list, V_list(:,4),'k.','LineWidth', 1);
legend('X Velocity', 'Y Velocity')
title('Velocity')


%% --- Scatter Plot: Step Size vs Distance to Sun ---
r_list = sqrt(V_var(:,1).^2 + V_var(:,2).^2);   % distance from planet to sun
h_step_list = diff(t_var);                      % step sizes from variable-step RK

% Note: r_list has one more element than h_step_list, so use r_list(1:end-1)
figure();
scatter(r_list(1:end-1), h_step_list, 25, 'filled');
xlabel('Distance to Sun, r');
ylabel('Step Size, h');
title('Variable Step: Step Size vs Distance to Sun');
grid on;