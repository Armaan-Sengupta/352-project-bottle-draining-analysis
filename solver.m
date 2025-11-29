clear; clc;

%% 1. Define Parameters
% Given values
rho_val = 998.21;       % Density (kg/m^3)
mu_val = 0.001003;      % Dynamic Viscosity mu (formerly N) [kg/m.s]
D_val = 0.11045;        % Tank Diameter (m)
d_val = 0.0065;         % Pipe Diameter (m)
K2_val = 0.78;          % Loss Coefficient
h0_val = 0.1365;        % Initial Height (m)

L_val = 0.5;            % Pipe Length (m) 
epsilon_val = 0.5e-6;   % Pipe Roughness (m) 
g_val = 9.81;           % Gravity (m/s^2)

%% 2. Symbolic Setup 
% Defining 'mu' specifically as requested
syms h(t) dh v2
syms rho mu D d K2 L epsilon g real

% Continuity Equation: v2 related to dh (rate of change of height)
% Note: dh is negative as tank drains, so we use -dh to get positive velocity
eq_continuity = v2 == (D/d)^2 * (-dh); 

% Reynolds Number
% Note: Standard pipe flow uses diameter 'd' for Re, not length 'L'
Re = (rho * v2 * d) / mu;

% Friction Factor (Haaland Equation approximation)
% Used for explicit symbolic representation instead of implicit Colebrook
f_sym = (-1.8 * log10( (epsilon/d)/3.7 + 6.9/Re ))^(-2);

% Energy Equation (Bernoulli with Head Loss)
% LHS: Total Energy at free surface (z=h, v~dh)
% RHS: Kinetic Energy at exit + Major Loss (Friction) + Minor Loss (K2)
LHS = g*h + 0.5*(dh)^2;
RHS = 0.5*v2^2 + g * f_sym * (v2^2 / (2*g)) * (L/d) + K2 * (v2^2 / (2*g));

% Assemble ODE
ODEq = LHS == RHS;

% Substitute v2 to get the final ODE in terms of h and dh
ODEq_final = subs(ODEq, v2, rhs(eq_continuity));

disp('Symbolic Non-Linear ODE formulation:');
pretty(ODEq_final);

%% 3. Numerical Solution
t_span = [0 60]; % Time duration

% Pack parameters
p.rho = rho_val; p.mu = mu_val; p.D = D_val; p.d = d_val;
p.K2 = K2_val; p.L = L_val; p.eps = epsilon_val; p.g = g_val;

% Solve using ode15s (better for stiff/differential-algebraic problems)
opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-6);
[t_sol, h_sol] = ode15s(@(t,h) tank_ode(t, h, p), t_span, h0_val, opts);

% Plotting
figure;
plot(t_sol, h_sol, 'b-', 'LineWidth', 2);
yline(0, 'r--', 'Empty Tank');
xlabel('Time (s)');
ylabel('Height (m)');
title('Tank Discharge: Height vs Time');
grid on;

%% 4. Helper Function for ODE Solver
function dhdt = tank_ode(~, h, p)
    % Physics Check: If height is zero or negative, flow stops
    if h <= 1e-4
        dhdt = 0;
        return;
    end
    
    % Define the implicit equation for velocity v2 at current height h
    % We solve for v2 first, then convert to dhdt
    % Energy Balance: g*h = (1 + f*L/d + K2) * v^2 / 2g
    % (Neglecting the small surface velocity head 0.5*dh^2 for stability in the guess)
    
    balance_func = @(v) solve_energy_balance(v, h, p);
    
    % Initial Guess for v2 (Torricelli with no losses)
    v_guess = sqrt(2 * p.g * h);
    
    % Solve for exit velocity v2
    try
        v_sol = fzero(balance_func, v_guess);
    catch
        v_sol = 0;
    end
    
    % Convert v2 back to dh/dt using Continuity
    % dh/dt = - v2 * (d/D)^2
    dhdt = -v_sol * (p.d / p.D)^2;
end

function err = solve_energy_balance(v, h, p)
    % Calculate Re and f based on current velocity guess v
    if v <= 0
        f = 0; 
    else
        Re = (p.rho * v * p.d) / p.mu;
        if Re < 2300
            f = 64/Re; 
        else
            term = (p.eps/p.d)/3.7 + 6.9/Re;
            f = (-1.8 * log10(term))^(-2);
        end
    end
    
    % Calculate Head Loss terms
    h_friction = f * (p.L / p.d) * (v^2 / (2*p.g));
    h_minor = p.K2 * (v^2 / (2*p.g));
    h_kinetic = v^2 / (2*p.g);
    
    % Error = Available Head - (Velocity Head + Losses)
    % Note: We neglect the surface kinetic head (0.5*dh^2) here as it 
    % creates a circular dependency and is mathematically negligible for D >> d
    err = h - (h_kinetic + h_friction + h_minor);
end