%this function computes the orbit of a planet about a sun
%the sun is assumed to located at the origin (0,0)
%and motion is restricted to the x-y plane
%INPUTS:
%t_list: a list of times to compute the position & velocity of the planet
%V0: initial condition.  V0 is a column vector consisting
%   of the initial position and velocity of the planet:
%       V0 = [x(0); y(0); dx/dt(0); dy/dt(0)]
%orbit_params: a struct describing the system parameters
%   orbit_params.m_sun: mass of the sun
%   orbit_params.m_planet: mass of the planet
%   orbit_params.G: gravitational constant
%       Force = -m_planet*m_sun*G/r^2
%OUTPUTS:
%V_list: the state of the planet at each time in t_list
%   if t_list is a list, then V_list is a Nx4 MATRIX
%   where each ROW has the form [x_i,y_i,dx/dt_i,dy/dt_i]
%   corresponding to the values at the ith time
%   if t_list is a SCALAR (i.e. t_list = t),
%   then V_list is a COLUMN VECTOR of the form:
%   [x(t); y(t); dx/dt(t); dy/dt(t)]
%NOTES:
%This function needs all the other functions in this file to run
%I HIGHLY RECOMMEND JUST SAVING THIS FUNCTION IN ITS OWN FILE
%DON'T COPY AND PASTE INTO YOUR CODE! IT'S NOT WORTH IT!
%
%USAGE EXAMPLE:
%At the bottom of this file is a function usage_example()
%which shows how to use compute_planetary_motion(...)
%You can start from there and then tweak it.
function V_list = compute_planetary_motion(t_list,V0,orbit_params)
    m_sun = orbit_params.m_sun;
    m_planet = orbit_params.m_planet;
    G = orbit_params.G;

    i = sqrt(-1);

    x0 = V0(1);
    y0 = V0(2);
    vx0 = V0(3);
    vy0 = V0(4);

    Q1 = abs(x0*vy0-y0*vx0);
    Q2 = abs(1/2*(vx0^2+vy0^2)-m_sun*G/sqrt(x0^2+y0^2));

    tau = Q1/Q2;
    l_scale = sqrt(Q1^2/Q2);
    velocity_scale = sqrt(Q2);

    x0=x0/l_scale;
    y0=y0/l_scale;
    vx0 = vx0/velocity_scale;
    vy0 = vy0/velocity_scale;

    t_list = t_list/tau;


    G = m_sun*G/(Q1*sqrt(Q2));
    m_sun = 1;
    m_planet = 1;


    H0 = m_planet*(x0*vy0-y0*vx0);

    dArea_dt = .5*H0/m_planet;


    r0 = abs(x0+y0*i);
    theta0 = angle(x0+y0*i);

    rdot0 = (vx0*x0+vy0*y0)/r0;
    theta_dot0 = H0/(m_planet*r0^2);

    alpha = m_sun*m_planet^2*G/H0^2;

    u0=1/r0;

    M = [-cos(theta0),sin(theta0);sin(theta0),cos(theta0)];
    B = [u0-alpha;-rdot0*u0^2/theta_dot0];
    Q = M\B;

    A_polar = abs(Q(1)+Q(2)*i);
    phi = angle(Q(1)+Q(2)*i);

    num_pts = length(t_list);
    V_list = zeros(num_pts,4);
    

    if 0<=alpha-A_polar && alpha-A_polar<1e-4
        d = 2/(alpha+A_polar);
        h0 = sin(theta0+phi)/(alpha+A_polar*cos(theta0+phi));
        area0 = kepler_func_parabola(h0,d,0);
    else
        a = alpha/abs(alpha^2-A_polar^2);
        b = sqrt(1/abs(alpha^2-A_polar^2));
        c = A_polar/abs(A_polar^2-alpha^2);

        if alpha>A_polar
            ellipse_area = pi*b*a;
            area0 = kepler_func_ellipse(theta0+phi,a,b,0);
        else
            h0 = (b^2/(a-c*cos(theta0+phi)))*sin(theta0+phi);
            area0 = kepler_func_hyperbola(h0,a,b,0);
        end
    end

    for n = 1:num_pts
        t = t_list(n);

        if 0<=alpha-A_polar && alpha-A_polar<1e-4
            target_area = -t*dArea_dt+area0;

            error_func = @(h_in) kepler_func_parabola(h_in,d,target_area);
            h_out = mini_secant_method(error_func,0);

            l_out = .5*h_out.^2/d -d/2;

            x_out = cos(-phi)*l_out-sin(-phi)*h_out;
            y_out = sin(-phi)*l_out+cos(-phi)*h_out;

            theta_out = angle(x_out+i*y_out);
            r_out = abs(x_out+i*y_out);
        else
            if alpha>A_polar
                target_area = t*dArea_dt+area0;
                target_area = mod(target_area,ellipse_area);
    
                if target_area>ellipse_area/2
                    target_area = target_area-ellipse_area;
                end
    
                error_func = @(theta_in) kepler_func_ellipse(theta_in,a,b,target_area);
                theta_out = mini_secant_method(error_func,0)-phi;
    
                r_out =  1/(alpha-A_polar*cos(theta_out+phi));

                x_out = r_out*cos(theta_out);
                y_out = r_out*sin(theta_out);
            else
                target_area = -t*dArea_dt+area0;

                error_func = @(h_in) kepler_func_hyperbola(h_in,a,b,target_area);
                h_out = mini_secant_method(error_func,0);
    
                l_out = a*sqrt(1+h_out^2/b^2)-c;
                
    
                x_out = cos(-phi)*l_out-sin(-phi)*h_out;
                y_out = sin(-phi)*l_out+cos(-phi)*h_out;
    
                theta_out = angle(x_out+i*y_out);
                r_out = abs(x_out+i*y_out);
            end
        end

        dtheta_dt_out = H0/(m_planet*r_out^2);
        drdt_out = -A_polar*sin(theta_out+phi)*dtheta_dt_out/(alpha-A_polar*cos(theta_out+phi))^2;

        vx_out = drdt_out*cos(theta_out)-dtheta_dt_out*r_out*sin(theta_out);
        vy_out = drdt_out*sin(theta_out)+dtheta_dt_out*r_out*cos(theta_out);

        V_list(n,:)=[x_out,y_out,vx_out,vy_out];
    end


    V_list = [V_list(:,1)*l_scale, V_list(:,2)*l_scale,...
            V_list(:,3)*velocity_scale, V_list(:,4)*velocity_scale];

    if isscalar(t_list)
        V_list = V_list';
    end

end

function error_val = kepler_func_parabola(y,d,target_area)
    error_val = y^3/(12*d)+y*d/4 - target_area;
end

function error_val = kepler_func_hyperbola(y,a,b,target_area)
    c = sqrt(a^2+b^2);
    q = y/b;
    hyperbolic_angle = log(q+sqrt(1+q^2));

    error_val = y*c/2-a*b*hyperbolic_angle/2 - target_area;
end

function error_val = kepler_func_ellipse(theta,a,b,target_area)
    i = sqrt(-1);

    c = sqrt(a^2-b^2);
    r = b^2/(a-c*cos(theta));

    x = (r*cos(theta)-c)/a;
    y = (r*sin(theta))/b;
    psi = angle(x+i*y);

    error_val = a*b*psi/2 + (c/2)*r*sin(theta) - target_area;
end

function x_out = mini_secant_method(func_in,xa)
    my_tol = 1e-13;

    xb = xa+1e-2;

    ya = func_in(xa);
    yb = func_in(xb);

    count = 0;
    while count<20 && abs(yb)>my_tol
        count = count+1;

        xc = (xa*yb-xb*ya)/(yb-ya);
        yc = func_in(xc);

        xa = xb;
        ya = yb;

        xb = xc;
        yb = yc;
    end

    x_out = xb;
end

%example for how to use compute_planetary_motion(...)
function usage_example()
    orbit_params = struct();
    orbit_params.m_sun = 1;
    orbit_params.m_planet = 1;
    orbit_params.G = 40;
    x0 = 8;
    y0 = 0;
    dxdt0 = 0;
    dydt0 = 1.5;
    
    V0 = [x0;y0;dxdt0;dydt0];
    t_range = linspace(0,30,100);
    V_list = compute_planetary_motion(t_range,V0,orbit_params);
    
    
    axis equal; axis square;
    axis([-20,20,-20,20])
    hold on
    plot(0,0,'ro','markerfacecolor','r','markersize',5);
    plot(V_list(:,1),V_list(:,2),'k');
end