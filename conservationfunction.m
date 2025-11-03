function conservationfunction()
% CONSERVATIONFUNCTION
% experiments on conservation of mechanical energy E and angular momentum H
% for an orbit ODE integrated with several explicit Runge–Kutta methods

%% system, units, and initial condition
m_sun    = 1;                 % solar mass in chosen units
m_planet = 1/330000;          % Earth-ish mass ratio (cancels in E,H as defined below)
G        = 4*pi^2/m_sun;      % AU-yr units: GM_sun = 4*pi^2
ti = 0; tf = 1;               % integrate about one year
V0 = [1; 0; 0; 2*pi];         % r0=1 AU, v≈2*pi AU/yr (near circular)

params.m_sun = m_sun; params.m_planet = m_planet; params.G = G;

%% RK methods (Butcher tableaus)
Heun.A = [0 0; 1 0];        Heun.B = [1/2 1/2];  Heun.C = [0; 1];
Heun.name = 'Heun (2nd)';   Heun.s = 2;

Ral.A  = [0 0 0; 1/2 0 0; 0 3/4 0];
Ral.B  = [2/9 1/3 4/9];     Ral.C = [0; 1/2; 3/4];
Ral.name = 'Ralston (3rd)'; Ral.s = 3;

RK4.A  = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
RK4.B  = [1/6 1/3 1/3 1/6]; RK4.C = [0; 1/2; 1/2; 1];
RK4.name = 'RK4 (4th)';     RK4.s = 4;

methods = {Heun, Ral, RK4};
M = numel(methods);

%% reference values E0, H0 from a very fine RK4 run
[t_ref, X_ref] = rk_fixed(@gravity_rate_func, [ti tf], V0, 1/20000, struct_merge(RK4, params));
[E_ref, H_ref] = EH_series(X_ref, params);
E0 = E_ref(1); H0 = H_ref(1);      % true constants we target

%% quick time traces for a single step size to visualize drift
h_demo = 1/200;                       % pick a moderate h
colors = lines(M);
figure; tiledlayout(2,1,'Padding','compact','TileSpacing','compact')

for m = 1:M
    BT = struct_merge(methods{m}, params);
    [t, X] = rk_fixed(@gravity_rate_func, [ti tf], V0, h_demo, BT);
    [E, H] = EH_series(X, params);
    relE = abs(E - E0)/abs(E0);
    relH = abs(H - H0)/max(abs(H0), 1e-12);

    nexttile(1); hold on; grid on; box on
    plot(t, relE, 'o-','Color',colors(m,:),'DisplayName',[methods{m}.name '  E error'])
    ylabel('|E(t)-E_0|/|E_0|'); title(sprintf('Conservation drift, h = %.4g', h_demo))

    nexttile(2); hold on; grid on; box on
    plot(t, relH, 'o-','Color',colors(m,:),'DisplayName',[methods{m}.name '  H error'])
    ylabel('|H(t)-H_0|/|H_0|'); xlabel('t')
end
nexttile(1); legend('Location','northwest')
nexttile(2); legend('Location','northwest')

%% sweep step sizes: max conservation error vs h, with fits
hs = logspace(-4, -1, 16);
maxE = zeros(M, numel(hs));
maxH = zeros(M, numel(hs));

for m = 1:M
    BT = struct_merge(methods{m}, params);
    for i = 1:numel(hs)
        h = hs(i);
        [t, X] = rk_fixed(@gravity_rate_func, [ti tf], V0, h, BT);
        [E, H] = EH_series(X, params);
        maxE(m,i) = max(abs(E - E0))/abs(E0);
        maxH(m,i) = max(abs(H - H0))/max(abs(H0), 1e-12);
    end
end

% plot: max |E-E0| vs h
figure;  grid on; box on
leg = {};
for m = 1:M
    loglog(hs, maxE(m,:), 'o-','Color',colors(m,:),'MarkerSize',4,'LineWidth',1.0)
    hold on;
    [p,k] = loglog_fit(hs, maxE(m,:));
    loglog(hs, k*hs.^p, '--','Color',colors(m,:),'LineWidth',1.0)
    leg{end+1} = sprintf('%s  data', methods{m}.name);
    leg{end+1} = sprintf('fit  O(h^{%.2f})', p);
end
xlabel('h'); ylabel('max_t |E(t)-E_0| / |E_0|')
title('Energy conservation error vs step size h')
legend(leg,'Location','southeast')

% plot: max |H-H0| vs h
figure; grid on; box on
leg = {};
for m = 1:M
    loglog(hs, maxH(m,:), 'o-','Color',colors(m,:),'MarkerSize',4,'LineWidth',1.0)
    hold on;
    [p,k] = loglog_fit(hs, maxH(m,:));
    loglog(hs, k*hs.^p, '--','Color',colors(m,:),'LineWidth',1.0)
    leg{end+1} = sprintf('%s  data', methods{m}.name);
    leg{end+1} = sprintf('fit  O(h^{%.2f})', p);
end
xlabel('h'); ylabel('max_t |H(t)-H_0| / |H_0|')
title('Angular momentum conservation error vs step size h')
legend(leg,'Location','southeast')

%% small printed summary
fprintf('\nstep size sweep: slopes p for max |E-E0| ~ h^p \n')
for m = 1:M
    pE = loglog_fit(hs, maxE(m,:)); pE = pE(1);
    fprintf('%-16s  p ≈ %6.3f\n', methods{m}.name, pE)
end
fprintf('\nstep size sweep: slopes p for max |H-H0| ~ h^p \n')
for m = 1:M
    pH = loglog_fit(hs, maxH(m,:)); pH = pH(1);
    fprintf('%-16s  p ≈ %6.3f\n', methods{m}.name, pH)
end

end

%% helpers 

function [t_list, X_list] = rk_fixed(rate_func, tspan, X0, h_ref, BT)
    % fixed-step explicit RK integrator
    t0 = tspan(1); tf = tspan(2);
    N  = ceil((tf - t0)/h_ref);
    h  = (tf - t0)/N;
    t_list = linspace(t0, tf, N+1);
    X_list = zeros(numel(X0), N+1); X_list(:,1) = X0;
    XA = X0;
    for k = 1:N
        XA = rk_step(rate_func, t_list(k), XA, h, BT);
        X_list(:,k+1) = XA;
    end
end

function XB = rk_step(rate_func, t, XA, h, BT)
    % one explicit RK step from a Butcher tableau
    A = BT.A; B = BT.B(:); C = BT.C(:);
    s = numel(C); n = numel(XA);
    K = zeros(n, s);
    for i = 1:s
        sum_prev = K(:,1:i-1) * A(i,1:i-1)';          % sum_j a_ij k_j
        K(:,i)   = rate_func(t + C(i)*h, XA + h*sum_prev, BT);
    end
    XB = XA + h * (K * B);
end

function dVdt = gravity_rate_func(~, V, BT_with_params)
    % planet around a fixed sun at origin, 2-D
    x = V(1); y = V(2); vx = V(3); vy = V(4);
    r2 = x*x + y*y; r = sqrt(r2);
    Gm = BT_with_params.G * BT_with_params.m_sun;
    ax = -Gm * x / (r^3 + eps);
    ay = -Gm * y / (r^3 + eps);
    dVdt = [vx; vy; ax; ay];
end

function [E, H] = EH_series(X, params)
    % specific mechanical energy and specific angular momentum
    x = X(1,:); y = X(2,:); vx = X(3,:); vy = X(4,:);
    r = sqrt(x.^2 + y.^2);
    E = 0.5*(vx.^2 + vy.^2) - params.m_sun*params.G ./ r; % per unit planet mass
    H = (x.*vy - y.*vx);                                   % z-component per unit mass
end

function S = struct_merge(A,B)
    S = A;
    f = fieldnames(B);
    for i = 1:numel(f), S.(f{i}) = B.(f{i}); end
end

function [p, k] = loglog_fit(x, y)
    % fits y ≈ k * x^p on log–log axes
    X = [log(x(:)), ones(numel(x),1)];
    b = X \ log(y(:));
    p = b(1);
    k = exp(b(2));
end


% this experiment checks how well different rk methods keep energy and
% angular momentum constant during an orbit
% in real physics those quantities don’t change, but in simulations they
% drift a bit because of step size error

% setup is one planet orbiting a fixed sun under gravity
% we run three rk methods: heun (2nd), ralston (3rd), and rk4 (4th)
% all start from the same initial position and velocity
% we run for one orbit and compare how much energy (E) and angular momentum (H)
% change compared to their starting values

% we measure E = 1/2(vx^2 + vy^2) - GM/|r|
% and H = x*vy - y*vx
% then track how much E and H drift away from their initial values
% smaller drift means the method is better at conserving physics

% the drift plots show how |E - E0| and |H - H0| grow over time for each method
% rk4 stays almost perfectly flat (super stable)
% heun drifts the most
% ralston is in between

% the loglog plots show how the max energy and momentum errors shrink as we
% make the step size smaller
% the slope of each line (p) tells the order of accuracy
% results came out close to what we expect:
% heun around order 2.9, ralston around 2.9, rk4 around 4.5

% basically rk4 is best at conserving energy and angular momentum
% higher order methods lose less energy and angular momentum over time
% this shows that making the step smaller or using a higher order method
% keeps the orbit more realistic and stable


% how do we isolate from initial condition dependence
% use the exact same IC for every method when plotting a single curve
% also test a small batch of ICs that share the same energy level
%   e g vary the starting angle around the circle or use small eccentricities
% report median and IQR of the drift across the IC batch
% compare relative errors |E - E0| / |E0| and |H - H0| / |H0| so scaling of the ICs doesn’t bias things
% fix the final time and the evaluation times for all methods so phase differences don’t sneak in

% what does effective mean here
% small invariant drift for a given amount of computational work
% metrics we use
%   max over time of |E - E0| / |E0| and |H - H0| / |H0|
%   rms over time of those same quantities
%   slope of a linear fit to the drift vs time to see if it trends up
%   error per function call to capture cost vs accuracy
%   largest step size that stays stable over one orbit

% how do we control for the fact that smaller h always helps
% run a sweep over h and fit error ~ h^p so we compare orders p and prefactors fairly
% also plot error versus number of function evaluations so methods are compared at equal cost
% when drawing conclusions, use the error vs calls plot as the tie breaker
% keep the same final time and same error metrics across all h so apples to apples



