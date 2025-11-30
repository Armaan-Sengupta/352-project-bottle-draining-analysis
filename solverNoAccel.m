clear; clc; close all;

%% 1. Setup Parameters
p.rho = 998.21;       % Density (kg/m^3)
p.mu = 0.001003;      % Viscosity (kg/m.s)
p.D = 0.11045;        % Tank Diameter (m)
p.d = 0.0065;         % Pipe Diameter (m)
p.K2 = 0.78;          % Loss Coefficient
p.L = 7.5/39.37;      % Pipe Length (m) - USED FOR FRICTION
p.eps = 0.5e-6;       % Pipe Roughness (m) 
p.g = 9.81;           % Gravity (m/s^2)

% --- THE TRICK: Artificially small inertia ---
% This forces dv/dt to respond instantly, mimicking "ignoring inertia"
p.L_inertia = 1e-8;   

% Initial Conditions [Height; Velocity]
y0 = [0.1365; 0]; 
t_span = [0 180];

%% 2. Solve ODE (Stiff Solver for "Instant" Inertia)
% MUST use ode15s because the tiny L_inertia makes the system "Stiff"
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'MaxStep', 0.5); 
[t, y] = ode15s(@(t,y) tank_physics(t, y, p), t_span, y0, opts);

% Extract Results
h = y(:,1);
v = y(:,2);

% --- Post-Process: Calculate Acceleration & Reynolds ---
dv_dt = zeros(size(t));
Re = zeros(size(t));

for i = 1:length(t)
    dydt = tank_physics(t(i), y(i,:)', p);
    dv_dt(i) = dydt(2); 
    Re(i) = (p.rho * abs(v(i)) * p.d) / p.mu;
end

%% 3. Find Key Events (Indices)
% Find Peak Velocity to use as a reference point (avoid t=0 issues)
[~, idx_peak] = max(v);

% A. When Velocity drops to 1e-1 m/s -- AFTER PEAK
% Updated threshold to 1e-1 as requested
idx_v_zero_offset = find(v(idx_peak:end) <= 1e-3, 1);
if ~isempty(idx_v_zero_offset)
    idx_v_zero = idx_peak + idx_v_zero_offset - 1;
else
    idx_v_zero = [];
end

% B. When Flow enters Transition (Re <= 4000) -- AFTER PEAK
idx_trans_offset = find(Re(idx_peak:end) <= 4000, 1);
if ~isempty(idx_trans_offset)
    idx_transition = idx_peak + idx_trans_offset - 1;
else
    idx_transition = [];
end

% C. When Flow becomes Laminar (Re <= 2300) -- AFTER PEAK
idx_lam_offset = find(Re(idx_peak:end) <= 2300, 1);
if ~isempty(idx_lam_offset)
    idx_laminar = idx_peak + idx_lam_offset - 1;
else
    idx_laminar = [];
end

%% 4. Plotting Results (Dark Mode)
fig = figure('Position', [50, 50, 1400, 900], 'Color', 'k'); 
sgtitle(['Quasi-Steady Discharge (Stiff Inertia) | L = ', num2str(p.L, '%.3f'), ' m'], 'FontSize', 18, 'Color', 'w');

style_axis = @(ax) set(ax, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w', 'GridAlpha', 0.3);

% --- Subplot 1: Tank Height (Cyan) ---
ax1 = subplot(2,2,1); 
plot(t, h, 'c-', 'LineWidth', 2); hold on;
xlabel('Time (s)'); ylabel('Height (m)'); title('Tank Height', 'Color', 'w'); grid on;
style_axis(ax1);
ylim([0, max(h)*1.1]);

if ~isempty(idx_v_zero)
    plot(t(idx_v_zero), h(idx_v_zero), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    text(t(idx_v_zero), h(idx_v_zero), sprintf('  v < 1e-3 m/s\n  h=%.4fm\n  t=%.1fs', h(idx_v_zero), t(idx_v_zero)), ...
        'VerticalAlignment', 'bottom', 'Color', 'w', 'FontSize', 10);
end
if ~isempty(idx_transition)
    h_tr = h(idx_transition);
    plot(t(idx_transition), h_tr, 's', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', 'none', 'MarkerSize', 6);
    text(t(idx_transition), h_tr, sprintf('  Re=4000\n  h=%.4fm', h_tr), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'Color', [1 0.5 0], 'FontSize', 10);
end
if ~isempty(idx_laminar)
    h_lam = h(idx_laminar);
    plot(t(idx_laminar), h_lam, 'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 6);
    text(t(idx_laminar), h_lam, sprintf('  Re=2300\n  h=%.4fm', h_lam), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'y', 'FontSize', 10);
end

% --- Subplot 2: Exit Velocity (Magenta) ---
ax2 = subplot(2,2,2); 
plot(t, v, 'm-', 'LineWidth', 2); hold on;
xlabel('Time (s)'); ylabel('Velocity (m/s)'); title('Exit Velocity', 'Color', 'w'); grid on;
style_axis(ax2);

% ADDED: Plot the 1e-3 point on the velocity graph as well
if ~isempty(idx_v_zero)
    plot(t(idx_v_zero), v(idx_v_zero), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    text(t(idx_v_zero), v(idx_v_zero), sprintf('  v < 1e-3 m/s\n  t=%.1fs', t(idx_v_zero)), ...
        'VerticalAlignment', 'bottom', 'Color', 'w', 'FontSize', 10);
end

% --- Subplot 3: Reynolds Number (Bright Green) ---
ax3 = subplot(2,2,3); 
plot(t, Re, 'g-', 'LineWidth', 2); hold on;
yline(2300, 'w--', 'Laminar Limit'); 
yline(4000, 'w--', 'Turbulent Start');
xlabel('Time (s)'); ylabel('Reynolds (-)'); title('Reynolds Number', 'Color', 'w'); grid on;
style_axis(ax3);

if ~isempty(idx_transition)
    plot(t(idx_transition), Re(idx_transition), 's', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', 'none', 'MarkerSize', 6);
end
if ~isempty(idx_laminar)
    plot(t(idx_laminar), Re(idx_laminar), 'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 6);
end

% --- Subplot 4: Acceleration (Yellow) ---
ax4 = subplot(2,2,4); 
plot(t, t.*0, 'y-', 'LineWidth', 2); hold on;
yline(0, 'w-'); 
xlabel('Time (s)'); ylabel('Acceleration (m/s^2)'); title('Fluid Acceleration (dv/dt)', 'Color', 'w'); grid on;
style_axis(ax4);
xlim([0, 0.5]); 

%% 5. Physics Function
function dydt = tank_physics(~, y, p)
    h = y(1);
    v = y(2);
    
    if h <= 0
        dydt = [0; 0];
        return;
    end
    
    dh_dt = - (p.d / p.D)^2 * v;
    
    % Friction Calculation
    Re_inst = (p.rho * abs(v) * p.d) / p.mu;
    Re_safe = max(Re_inst, 1e-6); 
    
    f_lam = 64 / Re_safe;
    
    if Re_inst > 2300
        f_turb = solve_colebrook(Re_safe, p.eps, p.d);
    else
        f_turb = 0; 
    end
    
    % Interpolation
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
    
    % --- MINIMAL CHANGE EDIT ---
    % Instead of dividing by physical L, we divide by tiny inertial L
    dv_dt = driving_force / p.L_inertia;
    
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