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
t_span = [0, 6];
V0 = [8; 0; 0; 1.5];
h_ref = 0.2;

% Create butcher tableau
BT_struct = struct();
BT_struct.A = [0];
BT_struct.B = [1];
BT_struct.C = [0];

% Calculate RK solution
[t_list, V_list, h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in, t_span, V0, h_ref, BT_struct);

% Calculate actual solution
solution_V_list = zeros(length(t_list), length(V0));
for i = 1:length(t_list)
    solution_V_list(i, :) = compute_planetary_motion(t_list(i), V0, orbit_params)';
end

figure(1)
colors = ['b', 'k', 'm', 'r.'];
% plot RK solution
for i = 1:length(V0)
    plot(t_list, V_list(:,i), colors(i));
    hold on;
end
% plot actual solution
for i = 1:length(V0)
    plot(t_list, solution_V_list(:,i), colors(i) + "--");
    hold on;
end

xlabel('Time');
ylabel('Solution');
grid on;





%%
function compute_planetary_motion_example()
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    
    V0 = [x0;y0;dxdt0;dydt0];
    t_range = linspace(0,30,100);
    V_list = compute_planetary_motion(t_range,V0,orbit_params);
    
    
    axis equal; axis square;
    axis([-20,20,-20,20])
    hold on
    plot(0,0,'ro','markerfacecolor','r','markersize',5);
    plot(V_list(:,1),V_list(:,2),'k');
end