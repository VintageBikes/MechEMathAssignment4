clear; clc; close all;

%% --- Parameters ---
orbit_params = struct();
orbit_params.m_sun = 1;
orbit_params.m_planet = 1;
orbit_params.G = 40;

t_span = [0 16];            % simulation time
V0 = [8; 0; 0; 1.5];        % initial state [x; y; vx; vy]

% Bogacki-Shampine Butcher Tableau (4-stage, 3/2 embedded)
Bogacki = struct();
Bogacki.C = [0, 1/2, 3/4, 1];
Bogacki.B = [2/9, 1/3, 4/9, 0; 
             7/24, 1/4, 1/3, 1/8]; 
Bogacki.A = [0,0,0,0; 1/2,0,0,0; 0,3/4,0,0; 2/9,1/3,4/9,0];

p = 1;                    % order for error scaling
error_desired = 1e-6;     % local truncation error
h_ref = 0.01;             % reference step size
h_fixed = 0.01;           % fixed step size

%% --- Define rate function ---
rate_func_in = @(t, V) gravity_rate_func(t, V, orbit_params);

%% --- Adaptive RK Integration ---
[t_ad, V_ad, h_avg_ad, num_evals_ad, num_failed_ad] = ...
    explicit_RK_variable_step_integration(rate_func_in, t_span, V0, h_ref, Bogacki, p, error_desired);

step_failure_rate = num_failed_ad / (length(t_ad) + num_failed_ad - 1);

% Compute global truncation error
solution_ad = zeros(length(t_ad), length(V0));
for i = 1:length(t_ad)
    solution_ad(i,:) = compute_planetary_motion(t_ad(i), V0, orbit_params)';
end
global_error_ad = max(vecnorm(V_ad - solution_ad, 2, 2));

fprintf('Adaptive RK:\n');
fprintf('Global error = %.3e\n', global_error_ad);
fprintf('Average step size = %.5f\n', h_avg_ad);
fprintf('Function evaluations = %d\n', num_evals_ad);
fprintf('Step failure rate = %.3f\n\n', step_failure_rate);

%% --- Fixed-step RK Integration ---
[t_fix, V_fix, num_evals_fix] = explicit_RK_fixed_step_integration(rate_func_in, t_span, V0, h_fixed, Bogacki);

solution_fix = zeros(length(t_fix), length(V0));
for i = 1:length(t_fix)
    solution_fix(i,:) = compute_planetary_motion(t_fix(i), V0, orbit_params)';
end
global_error_fix = max(vecnorm(V_fix - solution_fix, 2, 2));

fprintf('Fixed-step RK:\n');
fprintf('Global error = %.3e\n', global_error_fix);
fprintf('Average step size = %.5f\n', h_fixed);
fprintf('Function evaluations = %d\n\n', num_evals_fix);

%% --- Plots ---
% Orbit (x vs y)
figure;
plot(solution_ad(:,1), solution_ad(:,2),'b-','LineWidth',1.5); hold on;
plot(V_ad(:,1), V_ad(:,2),'ro'); axis equal;
xlabel('x'); ylabel('y'); legend('Exact','Adaptive'); title('Planetary Orbit');

% Step size vs distance
h_list = diff(t_ad);
r_list = sqrt(V_ad(1:end-1,1).^2 + V_ad(1:end-1,2).^2);
figure;
scatter(r_list, h_list, 20, 'filled');
xlabel('Distance from Sun'); ylabel('Step size h'); title('Adaptive step size vs distance');

% Positions vs time
figure;
plot(t_ad, V_ad(:,1), 'r-o', t_ad, V_ad(:,2), 'b-o');
xlabel('Time'); ylabel('Position'); legend('x','y'); title('Position vs Time (Adaptive)');

% Velocities vs time
figure;
plot(t_ad, V_ad(:,3), 'r-o', t_ad, V_ad(:,4), 'b-o');
xlabel('Time'); ylabel('Velocity'); legend('vx','vy'); title('Velocity vs Time (Adaptive)');