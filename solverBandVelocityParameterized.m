clear; clc; close all;

%% 1. Setup Parameters
p.rho = 998.21;       
p.mu = 0.001003;      
p.D = 0.11045;        
p.d = 0.0065;         
p.K2 = 7.8;           
p.eps = 0.5e-6;       
p.g = 9.81;           
y0 = [0.1365; 0];     
t_span = [0 300];     % Simulation duration

% Tube Lengths to Sweep (from empirical data)
% We use these to define the bounds of our "Band"
lengths_in = [7.5, 6.5, 5.5, 4.5, 3.5, 2.5, 1.5, 0.5];
lengths_m = lengths_in / 39.37;

%% 2. Run Parameter Sweep & Interpolate
% We need a common time vector to compare velocities at the exact same "t"
t_common = linspace(0, 300, 1000); 
v_matrix = zeros(length(t_common), length(lengths_m));

fprintf('Running simulation for %d tube lengths...\n', length(lengths_m));
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'MaxStep', 0.1); 

for i = 1:length(lengths_m)
    p.L = lengths_m(i);
    
    % Solve ODE
    [t, y] = ode45(@(t,y) tank_physics(t, y, p), t_span, y0, opts);
    v = y(:,2);
    
    % Interpolate result onto the common time vector
    % This allows us to compare "apples to apples" at every time step
    v_matrix(:, i) = interp1(t, v, t_common, 'pchip', 0);
end

%% 3. Calculate Statistics (Min/Max Envelope)
v_min = min(v_matrix, [], 2); % Min velocity at each time step
v_max = max(v_matrix, [], 2); % Max velocity at each time step
v_mean = mean(v_matrix, 2);   % Average velocity

%% 4. Plotting Results (Dark Mode)
figure('Color', 'k', 'Position', [100, 100, 1000, 600]); 
hold on; grid on; box on;

% Style the Axes
ax = gca;
ax.Color = 'k';
ax.XColor = 'w';
ax.YColor = 'w';
ax.GridColor = 'w';
ax.GridAlpha = 0.3;

% --- Plot 1: The Shaded Band (Range) ---
% We create a polygon defined by the Min and Max curves
t_fill = [t_common, fliplr(t_common)];
v_fill = [v_min', fliplr(v_max')];

fill(t_fill, v_fill, 'c', 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
    'DisplayName', 'Velocity Range (All Lengths)');

% --- Plot 2: Individual Traces (Faint) ---
% Show the specific paths inside the band
for i = 1:length(lengths_m)
    plot(t_common, v_matrix(:,i), 'Color', [0 1 1 0.15], 'LineWidth', 0.5, ...
        'HandleVisibility', 'off'); % Don't show in legend
end

% --- Plot 3: Mean & Boundaries ---
plot(t_common, v_mean, 'w-', 'LineWidth', 2, 'DisplayName', 'Mean Velocity');
plot(t_common, v_max, 'c--', 'LineWidth', 1, 'DisplayName', 'Max Envelope (Shortest Tube starts here)');
plot(t_common, v_min, 'b--', 'LineWidth', 1, 'DisplayName', 'Min Envelope (Longest Tube starts here)');

% Formatting
xlabel('Time (s)', 'FontSize', 12, 'Color', 'w');
ylabel('Exit Velocity (m/s)', 'FontSize', 12, 'Color', 'w');
title('Velocity Envelope: Variation across Tube Lengths', 'FontSize', 14, 'Color', 'w');
legend('Location', 'northeast', 'FontSize', 10, 'TextColor', 'w', 'EdgeColor', 'w', 'Color', 'k');

% Zoom in on the interesting part (before everyone is drained)
xlim([0, 200]);
ylim([0, max(v_max)*1.1]);

hold off;

%% 5. Physics Function
function dydt = tank_physics(~, y, p)
    h = y(1);
    v = y(2);
    
    % Decay logic for empty tank to ensure smooth interpolation at the tail
    if h <= 0
        dydt = [0; -v]; 
        return;
    end
    
    dh_dt = - (p.d / p.D)^2 * v;
    
    Re_inst = (p.rho * abs(v) * p.d) / p.mu;
    Re_safe = max(Re_inst, 1e-6); 
    
    f_lam = 64 / Re_safe;
    
    if Re_inst > 2300
        f_turb = solve_colebrook(Re_safe, p.eps, p.d);
    else
        f_turb = 0; 
    end
    
    if Re_inst < 2300
        f = f_lam;
    elseif Re_inst > 4000
        f = f_turb;
    else
        w = (Re_inst - 2300) / (4000 - 2300);
        f = (1 - w) * f_lam + w * f_turb;
    end
    
    loss_friction = p.g * f * (p.L/p.d) * (v^2 / (2*p.g));
    loss_minor    = p.g * p.K2 * (v^2 / (2*p.g));
    
    driving_force = p.g*h + 0.5*dh_dt^2 - 0.5*v^2 - loss_friction - loss_minor;
    dv_dt = driving_force / p.L;
    
    dydt = [dh_dt; dv_dt];
end

%% 6. Helper: Iterative Colebrook Solver
function f = solve_colebrook(Re, epsilon, d)
    term = (epsilon/d)/3.7 + 6.9/Re;
    f_guess = (-1.8 * log10(term))^(-2);
    
    x = 1 / sqrt(f_guess); 
    tol = 1e-6; 
    max_iter = 10;
    
    A = (epsilon/d) / 3.7;
    B = 2.51 / Re;
    
    for i = 1:max_iter
        arg = A + B * x;
        if arg <= 0; arg = 1e-9; end 
        
        G = x + 2 * log10(arg);
        dG = 1 + (2 / log(10)) * (B / arg);
        
        x_new = x - G / dG;
        if abs(x_new - x) < tol
            x = x_new;
            break;
        end
        x = x_new;
    end
    f = 1 / x^2;
end