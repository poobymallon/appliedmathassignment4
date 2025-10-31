function cooper1017jahz()
    % system & ICs 
    m_sun    = 1;
    m_planet = 1/330000;
    G        = 4*pi^2/m_sun;      % 1 AU, 1 yr units
    ti = 0; tf = 1;               % one year keeps scales tidy
    V0 = [1; 0; 0; 2*pi];         % circular-ish: r0=1, v=2π AU/yr

    orbit_params.m_sun = m_sun;
    orbit_params.m_planet = m_planet;
    orbit_params.G = G;

    % RK methods (Butcher
    Heun.A = [0 0; 1 0];
    Heun.B = [1/2 1/2];
    Heun.C = [0; 1];
    Heun.name = 'Heun (2nd)';   Heun.s = 2;

    Ral.A = [0 0 0; 1/2 0 0; 0 3/4 0];
    Ral.B = [2/9 1/3 4/9];
    Ral.C = [0; 1/2; 3/4];
    Ral.name = 'Ralston (3rd)'; Ral.s = 3;

    RK4.A = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
    RK4.B = [1/6 1/3 1/3 1/6];
    RK4.C = [0; 1/2; 1/2; 1];
    RK4.name = 'RK4 (4th)';     RK4.s = 4;

    methods = {Heun, Ral, RK4};

    % reference (trut)
    % high-accuracy RK4 with very small step (no secant file needed)
    [t_ref,X_ref] = explicit_RK_fixed_step_integration(@gravity_rate_func, [ti tf], V0, 1/40000, struct_merge(RK4, orbit_params));
    truth = @(tq) interp1(t_ref', X_ref', tq, 'pchip')';   % 4x1 at tq

    % sweep over h 
    hs = logspace(-4, -1, 18);   % 1e-4 .. 1e-1
    M  = numel(methods);

    err_h   = zeros(M, numel(hs));
    evals_h = zeros(M, numel(hs));

    for m = 1:M
        BT = struct_merge(methods{m}, orbit_params);
        for i = 1:numel(hs)
            h = hs(i);
            [t, X, h_avg, num_evals] = explicit_RK_fixed_step_integration(@gravity_rate_func, [ti tf], V0, h, BT);
            Xtrue_tf = truth(t(end));                         % align exactly at tf
            err_h(m,i) = norm(X(1:2,end) - Xtrue_tf(1:2));   % position error
            evals_h(m,i) = num_evals;
        end
    end

    % plots with fit-lines 
    colors = lines(M);

    % A) global error vs step size h (with fit lines)
    figure;  grid on; box on
    leg = cell(1,2*M);
    for m = 1:M
        loglog(hs, err_h(m,:), 'o-', 'Color', colors(m,:), 'MarkerSize', 4, 'LineWidth', 1.0);
        hold on;
        [p,k] = loglog_fit(hs, err_h(m,:));
        loglog(hs, k*hs.^p, '--', 'Color', colors(m,:), 'LineWidth', 1.0);
        leg{2*m-1} = sprintf('%s data', methods{m}.name);
        leg{2*m}   = sprintf('fit: O(h^{%.2f})', p);
        p_table_h(m,1) = p; %#ok<AGROW>
    end
    xlabel('h'); ylabel('global position error at t_f');
    title('Global truncation error vs step size (with fits)');
    legend(leg, 'Location','southwest');

    % B) global error vs # of f-evals (with fit lines)
    figure;  grid on; box on
    leg = cell(1,2*M);
    for m = 1:M
        x = evals_h(m,:);
        y = err_h(m,:);
        loglog(x, y, 'o-', 'Color', colors(m,:), 'MarkerSize', 4, 'LineWidth', 1.0);
        hold on;
        [p,k] = loglog_fit(x, y);
        loglog(x, k*x.^p, '--', 'Color', colors(m,:), 'LineWidth', 1.0);
        leg{2*m-1} = sprintf('%s data', methods{m}.name);
        leg{2*m}   = sprintf('fit: O(calls^{%.2f})', p);
        p_table_calls(m,1) = p; 
    end
    xlabel('# of f calls'); ylabel('global position error at t_f');
    title('Global truncation error vs number of evaluations (with fits)');
    legend(leg, 'Location','southwest');

    % print the two small tables 
    fprintf('\n  Fitted slopes p for error vs h  \n');
    for m = 1:M
        fprintf('  %-14s  p ≈ %6.3f\n', methods{m}.name, p_table_h(m,1));
    end

    fprintf('\n  Fitted slopes p for error vs #evals   \n');
    for m = 1:M
        fprintf('  %-14s  p ≈ %6.3f\n', methods{m}.name, p_table_calls(m,1));
    end
end

% RK helpers 
function [XB, num_evals] = explicit_RK_step(rate_func_in,t,XA,h,BT)
    A = BT.A; B = BT.B(:); C = BT.C(:);
    s = numel(C); n = numel(XA);
    K = zeros(n,s);
    for i = 1:s
        sum_prev = K(:,1:i-1) * A(i,1:i-1)';                 % Σ a_ij k_j
        K(:,i)   = rate_func_in(t + C(i)*h, XA + h*sum_prev, BT);
    end
    XB = XA + h * (K * B);
    num_evals = s;
end

function [t_list, X_list, h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in, tspan, X0, h_ref, BT)
    ti = tspan(1); tf = tspan(2);
    N  = ceil((tf - ti)/h_ref);
    h_avg = (tf - ti)/N;
    t_list = linspace(ti, tf, N+1);
    X_list = zeros(numel(X0), N+1); X_list(:,1) = X0;
    num_evals = 0; XA = X0;
    for k = 1:N
        [XB, adds] = explicit_RK_step(rate_func_in, t_list(k), XA, h_avg, BT);
        num_evals = num_evals + adds;
        X_list(:,k+1) = XB; XA = XB;
    end
end

% dynamics 
function dVdt = gravity_rate_func(~,V,params)
    % V = [x; y; vx; vy], point-mass sun at origin
    x = V(1); y = V(2); vx = V(3); vy = V(4);
    r2 = x*x + y*y; r = sqrt(r2);
    Gm = params.G * params.m_sun;
    ax = -Gm * x / (r^3 + eps);     % -GM r_vec / r^3
    ay = -Gm * y / (r^3 + eps);
    dVdt = [vx; vy; ax; ay];
end

% utilities 
function S = struct_merge(A,B)
    S = A;
    f = fieldnames(B);
    for i = 1:numel(f), S.(f{i}) = B.(f{i}); end
end

function [p, k] = loglog_fit(x, y)
    % y ≈ k x^p  =>  log(y) = p log(x) + log(k)
    X = [log(x(:)), ones(numel(x),1)];
    b = X \ log(y(:));
    p = b(1); k = exp(b(2));
end
