function coopervariable()
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1/330000;
    orbit_params.G = 4*pi^2/orbit_params.m_sun;
    ti = 0;
    tf = 10;
    V0 = [1; 0; 0; 6.28;];

    DormandPrince = struct();
    DormandPrince.C = [0, 1/5, 3/10, 4/5, 8/9, 1, 1];
    DormandPrince.B = [35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0;...
    5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40];
    DormandPrince.A = [0,0,0,0,0,0,0;
    1/5, 0, 0, 0,0,0,0;...
    3/40, 9/40, 0, 0, 0, 0,0;...
    44/45, -56/15, 32/9, 0, 0, 0,0;...
    19372/6561, -25360/2187, 64448/6561, -212/729, 0, 0,0;...
    9017/3168, -355/33, 46732/5247, 49/176, -5103/18656, 0,0;...
    35/384, 0, 500/1113, 125/192, -2187/6784, 11/84,0];
    dorp = 4;

    tspan = [ti,tf];
    grav_wrapper = @(t, V) gravity_rate_func(t, V, orbit_params);

    %test working
    des_err = 1e-8;
    [t_list,X_list,~, ~] = explicit_RK_variable_step_integration(grav_wrapper, tspan, V0, 1/36.5, DormandPrince, 5, des_err);
    tforV = linspace(ti,tf, 366);
    V_list = compute_planetary_motion(tforV, V0, orbit_params);

    figure
    subplot(2,1,1)
    plot(t_list, X_list(1,:), "r--", "DisplayName", "Dormand-Prince")
    legend()
    title('x position over time - runge-kutta')
    subplot(2,1,2)
    plot(tforV, V_list(:,1), 'b', "DisplayName", "Secant Method")
    legend()
    title('x position over time - secant')

    figure
    subplot(2,1,1)
    plot(X_list(1,:), X_list(2,:), "r--", "DisplayName", "Dormand-Prince")
    legend()
    title('y vs x position of earth - runge-kutta')
    subplot(2,1,2)
    plot(V_list(:,1), V_list(:,2), 'b', "DisplayName", "Secant Method")
    legend()
    title('y vs x position of earth - secant')


    figure
    hold on
    plot(X_list(1,:), X_list(2,:), "r^", "DisplayName", "Dormand-Prince")
    plot(V_list(:,1), V_list(:,2), 'b', "DisplayName", "Secant Method")
    legend()
    title('y vs x position of earth')

    h_ref_list = logspace(-5,-1, 30);
    num_evals_list = [];
    h_avg_list = [];
    tr_error_list = [];

    % %in class local truncation figure
    % n_samples = 60;
    % h_ref_list = logspace(-3, 1, n_samples)
    % 
    % abs_diff_list = zeros(1, n_samples);
    % 
    % for i = 1:length(h_ref_list)
    %     h_ref = h_ref_list(i);
    %     V_list = compute_planetary_motion(ti+h_ref,V0,orbit_params);
    % 
    %     [XB1, ~, ~] = RK_step_embedded(grav_wrapper, ti, V0, h_ref, DormandPrince);
    % 
    %     abs_diff_list(i) = norm(V_list-V0);
    %     tr_error_list1(i) = norm(XB1-V_list);
    %     tr_error_list2(i) = norm(XB1-V_list);
    % end
    % 
    % figure
    % loglog(h_ref_list, abs_diff_list)
    % loglog(h_ref_list, tr_error_list1)
    % loglog(h_ref_list, tr_error_list2)



end

function [XB1, XB2, num_evals] = RK_step_embedded(rate_func_in,t,XA,h,BT_struct)
    as = BT_struct.A;
    bs = BT_struct.B;
    cs = BT_struct.C;
    s = length(cs);
    m = length(XA);
    ks = zeros(m, s);
    for i = 1:s
        SUMak = ks*(as(i,:)');
        ki = rate_func_in((t+cs(i)*h), (XA+h*SUMak));
        ks(:, i) = ki;
    end
    SUMbk1 = ks*bs(1,:)';
    SUMbk2 = ks*bs(2,:)';
    XB1 = XA + h*SUMbk1;
    XB2 = XA + h*SUMbk2;
    num_evals = s;
end

function [XB, num_evals, h_next, redo] = explicit_RK_variable_step(rate_func_in,t,XA,h,BT_struct,p,error_desired)
    alpha = 1.5;

    redo = false;
    [XB1, XB2, num_evals] = RK_step_embedded(rate_func_in,t,XA,h,BT_struct);
    epsc = norm(XB1-XB2);
    temp = (error_desired/epsc).^(1/p);
    h_next = min(0.9*temp, alpha)*h;
    XB = XB1;
    if error_desired < epsc
        redo = true;
    end
end

function [t_list, X_list, h_avg, num_evals] = variable_step_integration(rate_func_in, step_func, tspan, X0, h_ref, BT_struct, p, error_desired)
    ti = tspan(1); tf = tspan(2);
    N  = ceil((tf - ti)/h_ref); h = (tf - ti)/N;
    t_list = linspace(ti, tf, N+1);
    nx = numel(X0); X_list = zeros(nx, N+1); X_list(:,1) = X0;
    num_evals = 0; XA = X0;
    h_rec = []; 
    for k = 1:N
        redo = true;
        while redo == true
            h_rec(end+1,:) = h;
            [XB, add_evals, h, redo] = step_func(rate_func_in, t_list(k), XA, h, BT_struct, p, error_desired);
            num_evals = num_evals + add_evals;
        end
        X_list(:,k+1) = XB; XA = XB;
    end
    h_avg = mean(h_rec);
end

function [t_list,X_list,h_avg, num_evals] = explicit_RK_variable_step_integration(rate_func_in,tspan,X0,h_ref,BT_struct, p, error_desired)
    [t_list,X_list,h_avg, num_evals] = variable_step_integration(rate_func_in, @explicit_RK_variable_step, tspan, X0, h_ref, BT_struct, p, error_desired);
end

function dVdt = gravity_rate_func(t,V,orbit_params)
    xp = V(1);
    yp = V(2);
    r = [xp, yp];
    vxp = V(3);
    vyp = V(4);
    theta = atan(yp/xp);
    ms = orbit_params.m_sun;
    mp = orbit_params.m_planet;
    G = orbit_params.G;
    d = norm(r);
    acc = -ms*G/d;
    A = zeros(length(V), length(V));
    A(1, 3) = 1;
    A(2, 4) = 1;
    A(3, 1) = acc;
    A(4, 2) = acc;
    dVdt = A*V;
end

function DVdt = testfunc(t, X, params)
    DVdt = -5*X;
end