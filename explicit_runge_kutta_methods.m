% first-order methods
forward_euler = struct();
forward_euler.A = [0];
forward_euler.B = [1];
forward_euler.C = [0];

% second-order methods
explicit_midpoint_method = struct();
explicit_midpoint_method.A = [0 0; 1/2 0];
explicit_midpoint_method.B = [0 1];
explicit_midpoint_method.C = [0; 0.5];


heuns_method = struct();
heuns_method.A = [0 0; 1 0];
heuns_method.B = [1/2 1/2];
heuns_method.C = [0; 1];


ralstons_method = struct();
ralstons_method.A = [0 0; 2/3 0];
ralstons_method.B = [1/4 3/4];
ralstons_method.C = [0; 2/3];

% third-order methods
kuttas_third_order_method = struct();
kuttas_third_order_method.A = [0 0 0; 1/2 0 0; -1 2 0];
kuttas_third_order_method.B = [1/6 2/3 1/6];
kuttas_third_order_method.C = [0; 1/2; 1];


heuns_third_order_method = struct();
heuns_third_order_method.A = [0 0 0; 1/3 0 0; 0 2/3 0];
heuns_third_order_method.B = [1/4 0 3/4];
heuns_third_order_method.C = [0; 1/3; 2/3];


ralstons_third_order_method = struct();
ralstons_third_order_method.A = [0 0 0; 1/2 0 0; 0 3/4 0];
ralstons_third_order_method.B = [2/9 1/3 4/9];
ralstons_third_order_method.C = [0; 1/2; 3/4];


van_der_houwens_wrays_third_order_method = struct();
van_der_houwens_wrays_third_order_method.A = [0 0 0; 8/15 0 0; 1/4 5/12 0];
van_der_houwens_wrays_third_order_method.B = [1/4 0 3/4];
van_der_houwens_wrays_third_order_method.C = [0; 8/15; 2/3];


ssprk3 = struct();
ssprk3.A = [0 0 0; 1 0 0; 1/4 1/4 0];
ssprk3.B = [1/6 1/6 2/3];
ssprk3.C = [0; 1; 1/2];

% fourth-order methods
classic_fourth_order_method = struct();
classic_fourth_order_method.A = [
    0 0 0 0
    1/2 0 0 0
    0 1/2 0 0
    0 0 1 0
];
classic_fourth_order_method.B = [1/6 1/3 1/3 1/6];
classic_fourth_order_method.C = [0; 1/2; 1/2; 1];


three_eighths_rule_fourth_order_method = struct();
three_eighths_rule_fourth_order_method.A = [0 0 0 0; 1/3 0 0 0; -1/3 1 0 0; 0 0 1 0; 1 -1 1 0];
three_eighths_rule_fourth_order_method.B = [1/8 3/8 3/8 1/8];
three_eighths_rule_fourth_order_method.C = [0; 1/3; 2/3; 1];


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

% fifth-order methods
nystroms_fifth_order_method = struct();
nystroms_fifth_order_method.A = [
    0 0 0 0 0 0;
    1/3 0 0 0 0 0;
    4/25 6/25 0 0 0 0;
    1/4 -3 15/4 0 0 0;
    2/27 10/9 -50/81 8/81 0 0;
    2/25 12/25 2/15 8/75 0 0;
];
nystroms_fifth_order_method.B = [23/192 0 125/192 0 -27/64 125/192];
nystroms_fifth_order_method.C = [0; 1/3; 2/5; 1; 2/3; 4/5];