function coopervarmodified()
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1/330000;
    orbit_params.G = 4*pi^2/orbit_params.m_sun;
    ti = 0;
    % tf = 10;
    tf = 100;
    % V0 = [1; 0; 0; 6.28;];
    V0 = [1.8; 0; 0; 6.28;];

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

    % h0 = 1/36.5;
    h0 = 1;
    %test working
    des_err = 1e-5;
    [t_list,X_list,h_avg, num_evals, fail_rate, h_rec] = explicit_RK_variable_step_integration(grav_wrapper, tspan, V0, h0, DormandPrince, 5, des_err);
    % [t_list,X_list] = ode45(grav_wrapper, tspan, V0);
    % X_list = X_list';
    fail_rate
    tforV = linspace(ti,tf, 366);
    V_list = compute_planetary_motion(tforV, V0, orbit_params);

    figure
    % subplot(2,1,1)
    hold on
    plot(tforV, V_list(:,1), 'b', "DisplayName", "Secant Method")
    % plot(t_list, X_list(1,:), "r--", "DisplayName", "Dormand-Prince")

    plot(t_list, X_list(1,:), "k-", "DisplayName", "Dormand-Prince",'LineWidth',1);
    plot(t_list, X_list(1,:), "ro",'markerfacecolor','r', "DisplayName", "Dormand-Prince",'markersize',4);
    
    % xlim([min(t_list),max(t_list)])
    % legend()
    % title('x position over time - runge-kutta')
    % subplot(2,1,2)
    % plot(tforV, V_list(:,1), 'b', "DisplayName", "Secant Method")
    xlim([min(t_list),max(t_list)])
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
    plot(X_list(1,:), X_list(2,:), "k", "DisplayName", "Dormand-Prince",'LineWidth',1)
    plot(X_list(1,:), X_list(2,:), "ro", "DisplayName", "Dormand-Prince",'MarkerFaceColor','r','MarkerSize',4)
    plot(V_list(:,1), V_list(:,2), 'b', "DisplayName", "Secant Method")
    plot(0,0, 'ko','MarkerFaceColor','k','MarkerSize',4)
    axis equal
    legend()
    title('y vs x position of earth')

    h_ref_list = logspace(-5,-1, 30);
    num_evals_list = [];
    h_avg_list = [];
    tr_error_list = [];

    % time vs h - is step size changing dynamically?
    figure
    hold on
    plot(t_list(1:end-1), h_rec(1:end),'k','LineWidth',1)
    plot(t_list(1:end-1), h_rec(1:end),'ro','MarkerFaceColor','r','MarkerSize',5)

    %for your chosen method, do the following... (variable runge kutta)
    n_samples = 60;
    h_ref_list = logspace(-3, 1, n_samples);

    abs_diff_list = zeros(1, n_samples);
    approx_diff_list = zeros(1, n_samples);
    V_real = gravity_rate_func(ti,V0,orbit_params);
    XB1s = zeros(4, n_samples);
    XB2s = zeros(4, n_samples);
    for i = 1:length(h_ref_list)
        h_ref = h_ref_list(i);
        V_list = compute_planetary_motion(ti+h_ref,V0,orbit_params);
        % V_next = gravity_rate_func(t,V,orbit_params)
        [XB1, XB2, ~] = RK_step_embedded(grav_wrapper, ti, V0, h_ref, DormandPrince);
        % XB1s(:, i) = XB1';
        % XB2s(:, i) = XB2';
        abs_diff_list(i) = norm(V_list-V0);
        approx_diff_list(i) = norm(XB1-XB2);
        tr_error_list1(i) = norm(XB1-V_list);
        tr_error_list2(i) = norm(XB2-V_list);
    end

    figure
    loglog(h_ref_list, abs_diff_list, DisplayName='f(t+h)-f(t)')
    hold on
    loglog(h_ref_list, approx_diff_list, DisplayName='XB2-XB1')
    loglog(h_ref_list, tr_error_list1, DisplayName=['XB1 ' char(949) '_{local}'])
    loglog(h_ref_list, tr_error_list2, DisplayName=['XB2 ' char(949) '_{local}'])
    xlabel('Step Size (h)')
    ylabel('Errors')
    title('Local Differences vs. Step Size')
    legend(Location='southeast')

    [p_xb1,k_xb1] = loglog_fit(h_ref_list, tr_error_list1)
    [p_xb2,k_xb2] = loglog_fit(h_ref_list, tr_error_list2)
    [p_xb12,k_xb12] = loglog_fit(h_ref_list, approx_diff_list)

    
    % Plot the local truncation errors of XB1 and XB2 as a function of their difference, |XB1 − XB2|
    figure
    loglog(approx_diff_list, tr_error_list1, DisplayName=['XB1 ' char(949) '_{local}'])
    hold on
    loglog(approx_diff_list, tr_error_list2, DisplayName=['XB2 ' char(949) '_{local}'])
    loglog(approx_diff_list, approx_diff_list, 'k--', DisplayName='(XB1-XB2)')
    xlabel('XB1-XB2 (Approximated Difference)')
    ylabel('Differences to Compare')
    title(['XB1 and XB2 ', char(949), '_{local}' ' vs (XB1-XB2)'])
    legend(Location='southeast')

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

function [t_list, X_list, h_avg, num_fails, num_evals, h_rec] = variable_step_integration(rate_func_in, step_func, tspan, X0, h_ref, BT_struct, p, error_desired)
    ti = tspan(1); tf = tspan(2);
    N  = ceil((tf - ti)/h_ref); h = (tf - ti)/N;
    nx = numel(X0); 
    t_list = 1; t_list(1) = ti;
    X_list = zeros(nx, 1); X_list(:,1) = X0;
    num_evals = 0; XA = X0;
    h_rec = []; 
    num_fails = 0;
    tnow = ti;
    while tnow <= tf
        redo = true; 
        while redo == true
            h_prev = h;
            h = min(h,tf-tnow+1e-15);
            [XB, add_evals, h, redo] = step_func(rate_func_in, tnow, XA, h, BT_struct, p, error_desired);
            if redo == true
                num_fails = num_fails+1;
            end
            num_evals = num_evals + add_evals;
        end
        h_rec(end+1,:) = h_prev;
        tnow = tnow+h_prev;
        t_list(end+1) = tnow;
        X_list(:,end+1) = XB; XA = XB;
    end
    h_avg = mean(h_rec);
end

function [t_list,X_list,h_avg, num_evals, fail_rate, h_rec] = explicit_RK_variable_step_integration(rate_func_in,tspan,X0,h_ref,BT_struct, p, error_desired)
    [t_list,X_list,h_avg, num_fails, num_evals, h_rec] = variable_step_integration(rate_func_in, @explicit_RK_variable_step, tspan, X0, h_ref, BT_struct, p, error_desired);
    fail_rate = num_fails/num_evals;
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
    acc = -ms*G/d^3;
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

% log-log linear fit helper
function [p,k] = loglog_fit(x_regression,y_regression,varargin)
    if size(x_regression,1)==1, x_regression = abs(x_regression)'; end
    if size(y_regression,1)==1, y_regression = abs(y_regression)'; end
    if nargin==3
        fp = varargin{1}; N = length(x_regression); idx = (1:N).';
        mask = true(N,1);
        if isfield(fp,'min_index'), mask = mask & idx>=fp.min_index; end
        if isfield(fp,'max_index'), mask = mask & idx<=fp.max_index; end
        if isfield(fp,'min_xval'),  mask = mask & x_regression>=fp.min_xval; end
        if isfield(fp,'max_xval'),  mask = mask & x_regression<=fp.max_xval; end
        if isfield(fp,'min_yval'),  mask = mask & y_regression>=fp.min_yval; end
        if isfield(fp,'max_yval'),  mask = mask & y_regression<=fp.max_yval; end
        x_regression = x_regression(mask); y_regression = y_regression(mask);
    end
    Y = log(y_regression); X = [log(x_regression), ones(length(x_regression),1)];
    coeff = regress(Y, X); p = coeff(1); k = exp(coeff(2));
end
