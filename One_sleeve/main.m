clear all
close all
set(0, 'DefaultFigureRenderer', 'painters');
addpath("./functions");


sp.B = 2.8;        % Rod stiffness
sp.alpha = pi/2;   % Sleeve angle measured from the horizontal line
sp.m = 0.0;        % Mass at free tip
sp.gamma = 0.312;  % Distributed linear mass density
sp.L = 2;          % Rod total length
sp.ze = 0.0;       % Viscous dampling coefficient applied to mass m
sp.mu = 0;         % Friction coefficient
sp.g = [0; -9.8];  % Gravity vector

sp.F1 = @(t) [8*sin(4*pi*t); zeros(1,length(t))]; % Time-dependent Force applied to the free tip

sp.l0 = 1;    % Initial external length
sp.l0_t = 0;  % Initial rate of change of external length (-v1)

sp.tmax = 1;  % Max time for simulation
sp.dt = 5e-4; % timestep size
sp.Nel = 32;  % Number of elements

sp.max_iteration = 100; % Max number of iterations for the Newton-Raphson solver

sp.newmark_b1 = 0.255; % Newmark coefficient beta1
sp.newmark_b2 = 0.505; % Newmark coefficient beta2

sp.p = @(t) [0; 0];    % Sliding sleeve origin
sp.p_t = @(t) [0; 0];  % Sliding sleeve origin time derivative
sp.p_tt = @(t) [0; 0]; % Sliding sleeve origin second time derivative
sp.b = @(t) [ cos(sp.alpha); sin(sp.alpha) ]; % Sliding sleeve parallel vector
sp.b_t = @(t) [0; 0];
sp.b_tt = @(t) [0; 0];
sp.n = @(t) [ cos(sp.alpha + pi/2); sin(sp.alpha + pi/2) ]; % Sliding sleeve normal vector
sp.n_t = @(t) [0; 0];
sp.n_tt = @(t) [0; 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                          SOLVER SETUP                        %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 0;
t_prev = 0;

sp.map = setup_maps(sp.Nel);
sp.nodes = linspace(0,1,sp.Nel+1);

[U, Up, Upp] = initialize_dofs(t, sp);
Uinit = U;
U_prev = U;

s1 = sp.L - sp.l0;
v1 = 0;
a1 = 0;
s1_prev = s1;
v1_prev = v1;
a1_prev = a1;

fig = figure;
hold on
plot([U(sp.map.X(:,1)); U(sp.map.X(end,5))], [U(sp.map.X(:,2)); U(sp.map.X(end,6))], '-d' )
xlim([-0.05, 0.3])
ylim([-0.05, 0.3])
axis equal



Up_prev = Up;
Upp_prev = Upp;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                             SOLVER                           %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cnt = 0;
while( t<=sp.tmax )

    cnt = cnt + 1;
    dt = sp.dt;
    t_prev = t;
    t = t + dt;

    for kk=1:sp.max_iteration
        U_old = U;

        %%% Newton-Raphson iterations
        [UM, F, J11, Jxx2] = oss_newmark_matrices( U, U_prev, Up_prev, Upp_prev, s1, s1_prev, v1_prev, a1_prev, t, t_prev, sp );
        J = zeros(sp.map.Ndof+1);
        J(1:end-1, 1:end-1) = J11;
        J(end, 1:end-1) = jacobian_numerical( @(U_var) interface_condition(U_var, U_prev, Up_prev, Upp_prev, s1, s1_prev, v1_prev, a1_prev, t, t_prev, sp), U, 1e-4);
        J(:, end) = (nl_system_extended( [U; s1+1e-4], U_prev, Up_prev, Upp_prev, s1_prev, v1_prev, a1_prev, t, t_prev, sp ) - nl_system_extended( [U; s1-1e-4], U_prev, Up_prev, Upp_prev, s1_prev, v1_prev, a1_prev, t, t_prev, sp ) ) / (2*1e-4);

        FF = nl_system_extended( [U; s1], U_prev, Up_prev, Upp_prev, s1_prev, v1_prev, a1_prev, t, t_prev, sp );
        
        Jxx = J11(1:(4*(sp.Nel+1)), 1:(4*(sp.Nel+1)));
        Jxc = J11(1:(4*(sp.Nel+1)), 4*(sp.Nel+1) + [1:(sp.Nel+1+3)]);
        Jcx = J11(4*(sp.Nel+1) + [1:(sp.Nel+1+3)], 1:4*(sp.Nel+1));
        Jcc = J11(4*(sp.Nel+1) + [1:sp.Nel+1+3], 4*(sp.Nel+1) + [1:sp.Nel+1+3]);
        Jpp = J(4*(sp.Nel+1)+sp.Nel+1+3+[1:1], 4*(sp.Nel+1)+sp.Nel+1+3+[1:1]);
        J1p = J(1:(5*(sp.Nel+1)+3), 4*(sp.Nel+1)+sp.Nel+1+3+[1:1]);
        Jp1 = J(4*(sp.Nel+1)+sp.Nel+1+3+[1:1], 1:(5*(sp.Nel+1)+3));
        
        [Pe,Re,Ce] = equilibrate(Jxx);
        B = Re*Pe*Jxx*Ce;
        Jxxinv = Ce*( B\(Re*Pe) );

        J11sc  = Jcc - Jcx*Jxxinv*Jxc;
        J11sci = inv(J11sc);
        J11inv = [ Jxxinv + Jxxinv*Jxc*J11sci*Jcx*Jxxinv,  -Jxxinv*Jxc*J11sci;  -J11sci*Jcx*Jxxinv,  J11sci ];
        
        Jsc = Jpp - Jp1*J11inv*J1p;
        Jsci = inv(Jsc);
        Jinv11 = J11inv + J11inv*J1p*Jsci*Jp1*J11inv;
        Jinv12 = -J11inv*J1p*Jsci;
        Jinv21 = -Jsci*Jp1*J11inv;
        Jinv22 = Jsci;
        
        Jinv = [Jinv11, Jinv12; Jinv21, Jinv22];
        
        U = U - ( Jinv11*FF(1:end-1) + Jinv12*FF(end) );
        s1 = s1 - ( Jinv21*FF(1:end-1) + Jinv22*FF(end) );
        
        residual = norm(UM * U - F);
        step_size = norm(U-U_old);
        step_dir = (U-U_old); step_dir = step_dir/max(abs(step_dir));
        disp([num2str(kk),': res=', num2str(residual),' ,  step=', num2str(step_size)]);
        if residual < 1e-7*sp.map.Ndof && step_size < 1e-7*sp.map.Ndof
            break;
        end
    end
    
    Upp = 1/(dt*dt*sp.newmark_b1)*(U-U_prev) - 1/(dt*sp.newmark_b1)*Up_prev - (1-2*sp.newmark_b1)/(2*sp.newmark_b1)*Upp_prev;
    Up = Up_prev + (1-sp.newmark_b2)*dt*Upp_prev + sp.newmark_b2*dt*Upp;
    
    a1 = 1/(dt*dt*sp.newmark_b1)*(s1-s1_prev) - 1/(dt*sp.newmark_b1)*v1_prev - (1-2*sp.newmark_b1)/(2*sp.newmark_b1)*a1_prev;
    v1 = v1_prev + (1-sp.newmark_b2)*dt*a1_prev + sp.newmark_b2*dt*a1;
    
    tip_x(cnt) = U(sp.map.X(end,5));
    tip_y(cnt) = U(sp.map.X(end,6));
    tip_xv(cnt) = Up(sp.map.X(end,5));
    tip_yv(cnt) = Up(sp.map.X(end,6));
    time(cnt) = t;
    DOFS(:, cnt) = [U; s1];
    DOFS_p(:,cnt) = [Up; v1];
    DOFS_pp(:,cnt) = [Upp; a1];

    if mod(cnt,20)==0  %% Plotting
        clf(fig);
        direc = sp.b(0);
        subplot( 1,2, 1 )
        hold on
        plot([Uinit(sp.map.X(:,1)); Uinit(sp.map.X(end,5))], [Uinit(sp.map.X(:,2)); Uinit(sp.map.X(end,6))], '-*r' )
        plot([U(sp.map.X(:,1)); U(sp.map.X(end,5))], [U(sp.map.X(:,2)); U(sp.map.X(end,6))], '-*' )
        quiver( U(sp.map.X(:,1)), U(sp.map.X(:,2)), U(sp.map.X(:,3)), U(sp.map.X(:,4)) )
        plot(tip_x, tip_y, '-k');
        % xlim([-0.05, 1])
        % ylim([-0.05, 1])
        axis equal

        subplot( 1,2, 2 )
        hold on

        p(1)=plot(time, DOFS(end,:), 'g')
        p(2)=plot(time, tip_x-tip_x(1),'-r')
        p(3)=plot(time, tip_y-tip_y(1), '-b')

        xlabel('t')
        legend(p, {'s1', 'x(L)', 'y(L)'});
        drawnow
    end

    if s1<0 || s1>sp.L
        break;
    end

    U_prev = U;
    Up_prev = Up;
    Upp_prev = Upp;

    s1_prev = s1;
    v1_prev = v1;
    a1_prev = a1;
end


function Fval = nl_system_extended( U_ext, U_prev, Up_prev, Upp_prev, s1_prev, v1_prev, a1_prev, t, t_prev, sp )
    s1 = U_ext(end);
    U = U_ext(1:end-1);
    [UM, F] = oss_newmark_matrices( U, U_prev, Up_prev, Upp_prev, s1, s1_prev, v1_prev, a1_prev, t, t_prev, sp );
    IC = interface_condition(U, U_prev, Up_prev, Upp_prev, s1, s1_prev, v1_prev, a1_prev, t, t_prev, sp);
    Fval = [UM*U - F; IC];
end