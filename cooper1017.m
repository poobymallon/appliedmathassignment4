function cooper1017()
    m_sun = 1;
    m_planet = 1/330000;
    G = 4*pi^2/m_sun;
    ti = 0;
    tf = 10;
    V0 = [1; 0; 0; 6.28;];

    %Heun's method (s = 2)
    Heun_struct.A = [0, 0; 1, 0;];
    Heun_struct.B = [0.5, 0.5;];
    Heun_struct.C = [0; 1;];
    Heun_struct.m_sun = m_sun;
    Heun_struct.m_planet = m_planet;
    Heun_struct.G = G;


    %Ralston's third order method (s = 3)
    Ral_struct.A = [0,0,0;0.5,0,0;0,0.75,0;];
    Ral_struct.B = [2/9;1/3;4/9;];
    Ral_struct.C = [0;0.5;0.75;];
    Ral_struct.m_sun = m_sun;
    Ral_struct.m_planet = m_planet;
    Ral_struct.G = G;

    %og Runge-Kutta (s = 4)
    og_struct.A = [0,0,0,0;1/2,0,0,0;0,1/2,0,0;0,0,1,0;];
    og_struct.B = [1/6,1/3,1/3,1/6;];
    og_struct.C = [0;1/2;1/2;1;];
    og_struct.m_sun = m_sun;
    og_struct.m_planet = m_planet;
    og_struct.G = G;

    tspan = [ti,tf];
    [t_list,X_list,~, ~] = explicit_RK_fixed_step_integration(@gravity_rate_func,tspan,V0,1/365,og_struct);
    tforV = linspace(ti,tf, 365);
    V_list = compute_planetary_motion(tforV, V0, og_struct);

    figure
    subplot(2,1,1)
    plot(t_list, X_list(1,:), "r--", "DisplayName", "Runge-Kutta")
    legend()
    title('x position over time - runge-kutta')
    subplot(2,1,2)
    plot(tforV, V_list(:,1), 'b', "DisplayName", "Secant Method")
    legend()
    title('x position over time - secant')

    figure
    subplot(2,1,1)
    plot(X_list(1,:), X_list(2,:), "r--", "DisplayName", "Runge-Kutta")
    legend()
    title('y vs x position of earth - runge-kutta')
    subplot(2,1,2)
    plot(V_list(:,1), V_list(:,2), 'b', "DisplayName", "Secant Method")
    legend()
    title('y vs x position of earth - secant')


    figure
    hold on
    plot(X_list(1,:), X_list(2,:), "r--", "DisplayName", "Runge-Kutta")
    plot(V_list(:,1), V_list(:,2), 'b', "DisplayName", "Secant Method")
    legend()
    title('y vs x position of earth')


end

function [XB, num_evals] = explicit_RK_step(rate_func_in,t,XA,h,BT_struct)
    as = BT_struct.A;
    bs = BT_struct.B;
    cs = BT_struct.C;
    s = length(cs);
    m = length(XA);
    ks = zeros(m, s);
    for i = 1:s
        SUMak = ks*(as(i,:)');
        ki = rate_func_in((t+cs(i)*h), (XA+h*SUMak), BT_struct);
        ks(:, i) = ki;
    end
    SUMbk = ks*bs';
    XB = XA + h*SUMbk;
    num_evals = s;
end

function [t_list, X_list, h_avg, num_evals] = fixed_step_integration(rate_func_in, step_func, tspan, X0, h_ref, BT_struct)
    ti = tspan(1); tf = tspan(2);
    N  = ceil((tf - ti)/h_ref); h_avg = (tf - ti)/N;
    t_list = linspace(ti, tf, N+1);
    nx = numel(X0); X_list = zeros(nx, N+1); X_list(:,1) = X0;
    num_evals = 0; XA = X0;
    for k = 1:N
        [XB, add_evals] = step_func(rate_func_in, t_list(k), XA, h_avg, BT_struct);
        num_evals = num_evals + add_evals;
        X_list(:,k+1) = XB; XA = XB;
    end
end

function [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration(rate_func_in,tspan,X0,h_ref,BT_struct)
    [t_list,X_list,h_avg, num_evals] = fixed_step_integration(rate_func_in,@explicit_RK_step,tspan,X0,h_ref,BT_struct);
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