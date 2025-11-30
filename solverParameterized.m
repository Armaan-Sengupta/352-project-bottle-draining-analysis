clear; clc; close all;

%% 1. Input Empirical Data
% Data transcribed from the provided image
% Lengths in inches, Times in min:sec
emp_data = {
    7.5, '1:57';
    6.5, '1:52';
    5.5, '1:47';
    4.5, '1:40';
    3.5, '1:38';
    2.5, '1:34';
    1.5, '1:28';
    0.5, '1:25';
};

% Convert Empirical Data to Numeric Arrays (Meters and Seconds)
emp_L_in = [emp_data{:,1}]';
emp_time_str = emp_data(:,2);
emp_time_s = zeros(size(emp_L_in));

for i = 1:length(emp_time_str)
    parts = sscanf(emp_time_str{i}, '%d:%d');
    emp_time_s(i) = parts(1)*60 + parts(2);
end

% Convert length to meters for plotting consistency later
emp_L_m = emp_L_in / 39.37;

%% 2. Setup Simulation Parameters
p.rho = 998.21;       
p.mu = 0.001003;      
p.D = 0.11045;        
p.d = 0.0065;         
p.K2 = 0.0;           
p.eps = 0.5e-6;       
p.g = 9.81;           
y0 = [0.1365; 0];     
t_span = [0 300];     % Extended time to ensure draining for all lengths

% Simulate ONLY the empirical lengths
sim_L_in = emp_L_in; 
sim_drain_times = zeros(size(sim_L_in));

%% 3. Run Parameter Sweep
fprintf('Running simulation for %d specific tube lengths...\n', length(sim_L_in));

opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'MaxStep', 0.1); 

for i = 1:length(sim_L_in)
    % Set current length (convert inches to meters)
    p.L = sim_L_in(i) / 39.37;
    
    % Solve ODE
    [t, y] = ode45(@(t,y) tank_physics(t, y, p), t_span, y0, opts);
    h = y(:,1);
    v = y(:,2);
    
    % --- Find Drain Time (Strictly Velocity Based) ---
    [~, idx_peak] = max(v);
    
    % Check: Velocity drops below 1e-3 AFTER the peak
    idx_v_stop = find(v(idx_peak:end) <= 1e-3, 1);
    
    if ~isempty(idx_v_stop)
        idx_final = idx_peak + idx_v_stop - 1;
        sim_drain_times(i) = t(idx_final);
    else
        % Fallback if it never drains in time span
        sim_drain_times(i) = NaN; 
    end
end

%% 4. Linear Regressions (Simulated & Empirical)
% Filter out NaNs for Simulation
valid_idx = ~isnan(sim_drain_times);
sim_L_valid = sim_L_in(valid_idx);
sim_times_valid = sim_drain_times(valid_idx);

% Regression 1: Simulation
if length(sim_L_valid) > 1
    coeffs_sim = polyfit(sim_L_valid, sim_times_valid, 1);
    fitted_sim = polyval(coeffs_sim, sim_L_in);
    eqn_sim = sprintf('Sim Fit: T = %.2f L + %.2f', coeffs_sim(1), coeffs_sim(2));
else
    fitted_sim = zeros(size(sim_L_in));
    eqn_sim = 'Sim Fit: N/A';
end

% Regression 2: Empirical
coeffs_emp = polyfit(emp_L_in, emp_time_s, 1); 
fitted_emp = polyval(coeffs_emp, emp_L_in);
eqn_emp = sprintf('Emp Fit: T = %.2f L + %.2f', coeffs_emp(1), coeffs_emp(2));


%% 5. Plotting Results (Dark Mode)
figure('Color', 'k', 'Position', [100, 100, 1000, 600]); % Black Background
hold on; grid on; box on;

% Style the Axes
ax = gca;
ax.Color = 'k';
ax.XColor = 'w';
ax.YColor = 'w';
ax.GridColor = 'w';
ax.GridAlpha = 0.3;

% --- Plot Points ---
% Simulation Points (Cyan Circles)
plot(sim_L_in, sim_drain_times, 'co', 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'w', 'MarkerSize', 8, 'DisplayName', 'Sim Points');

% Empirical Points (Red Circles)
plot(emp_L_in, emp_time_s, 'ro', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'w', 'MarkerSize', 8, 'DisplayName', 'Emp Points');

% --- Plot Regression Lines ---
% Simulation Trend (Cyan Dashed)
plot(sim_L_in, fitted_sim, 'c--', 'LineWidth', 1.5, 'DisplayName', eqn_sim);

% Empirical Trend (Red Dashed)
plot(emp_L_in, fitted_emp, 'r--', 'LineWidth', 1.5, 'DisplayName', eqn_emp);

% Formatting
xlabel('Tube Length (inches)', 'FontSize', 12, 'Color', 'w');
ylabel('Time to Drain (s)', 'FontSize', 12, 'Color', 'w');
title('Drain Time Analysis: Simulation vs Experiment', 'FontSize', 14, 'Color', 'w');
legend('Location', 'northwest', 'FontSize', 10, 'TextColor', 'w', 'EdgeColor', 'w', 'Color', 'k');
xlim([0, 9]);

% Calculate Errors
if any(valid_idx)
    errors = sim_drain_times - emp_time_s'; 
    avg_abs_error = mean(abs(errors(~isnan(errors))));
    text(5, 50, sprintf('Mean Abs Error: %.2f s', avg_abs_error), 'FontSize', 10, 'Color', 'y', 'BackgroundColor', 'k', 'EdgeColor', 'y');
end

hold off;

%% 6. Physics Function
function dydt = tank_physics(~, y, p)
    h = y(1);
    v = y(2);
    
    % MODIFIED: If tank is empty, allow velocity to decay naturally
    % rather than freezing it. This allows v < 1e-3 check to work.
    if h <= 0
        dydt = [0; -v]; % Exponential decay of residual velocity
        return;
    end
    
    dh_dt = - (p.d / p.D)^2 * v;
    
    % Friction Calculation
    Re_inst = (p.rho * abs(v) * p.d) / p.mu;
    Re_safe = max(Re_inst, 1e-6); 
    
    f_lam = 64 / Re_safe;
    
    % Colebrook (Turbulent)
    if Re_inst > 2300
        f_turb = solve_colebrook(Re_safe, p.eps, p.d);
    else
        f_turb = 0; 
    end
    
    % Interpolation (Transition Region)
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
    
    % Bernoulli / Momentum Equation
    driving_force = p.g*h + 0.5*dh_dt^2 - 0.5*v^2 - loss_friction - loss_minor;
    dv_dt = driving_force / p.L;
    
    dydt = [dh_dt; dv_dt];
end

%% 7. Helper: Iterative Colebrook Solver
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