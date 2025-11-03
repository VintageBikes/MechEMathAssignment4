clear;
clc;
clf;
close all;

t = 0;
XA = [1; 2];
h = 0.1;

% Create Bogacki butcher tableau
Bogacki = struct();
Bogacki.C = [0,1/2, 3/4, 1];
Bogacki.B = [2/9, 1/3, 4/9, 0; 7/24, 1/4, 1/3, 1/8];
Bogacki.A = [0,0,0,0; 1/2,0,0,0; 0,3/4,0,0; 2/9,1/3, 4/9, 0];

p = 1;

error_desired = 1e-10;

t_span = [0, 7];

%[XB1, XB2, num_evals] = RK_step_embedded(@rate_func01,t,XA,h,Bogacki);

%[XB, num_evals2, h_next, redo] = explicit_RK_variable_step(@rate_func01,t,XA,h,Bogacki,p,error_desired);

[t_list,X_list,h_avg, num_evals] = explicit_RK_variable_step_integration(@rate_func01, t_span, XA, h, Bogacki, p, error_desired);

plot(t_list, X_list(:,1), '-k')
hold on;
plot(t_list, X_list(:,2))

function dXdt = rate_func01(t,X)
    dXdt = -5*X + 5*cos(t) - sin(t);
end

function X = solution01(t)
    X = cos(t);
end