clear all
close all
set(0, 'DefaultFigureRenderer', 'painters');
addpath("./functions");

sp.B = 0.15;       % Rod stiffness
sp.alpha = 0;      % Sleeve angle measured from the horizontal line
sp.gamma = 0.4;    % Distributed linear mass density
sp.L = 3;          % Rod total length
sp.c_distr = 0;    % Distributed drag coefficient
sp.mu = 0.0;       % Friction coefficient
sp.g = [0; -9.8];  % Gravity vector


sp.s10 = 1;  % Initial value for s1
sp.v10 = 0;  % Initial value for s1 time derivative
sp.s20 = 2;  % Initial value for s2
sp.v20 = 0;  % Initial value for s2 time derivative


sp.tmax = 15;  % Max time for simulation
sp.dt = 1e-3;  % timestep size
sp.Nel = 32;   % Number of elements

sp.max_iteration = 100; % Max number of iterations for the Newton-Raphson solver

sp.newmark_b1 = 0.255; % Newmark coefficient beta1
sp.newmark_b2 = 0.505; % Newmark coefficient beta2


sp.p1 = @(t) [0; 0]; % First sliding sleeve origin
sp.b1 = @(t) [ cos(sp.alpha); sin(sp.alpha) ]; % First sliding sleeve parallel vector
sp.n1 = @(t) [ cos(sp.alpha + pi/2); sin(sp.alpha + pi/2) ]; % First sliding sleeve normal vector

sp.p2 = @(t) [1; 0]; % Second sliding sleeve origin
sp.b2 = @(t) [ cos(sp.alpha); sin(sp.alpha) ]; % Second sliding sleeve parallel vector
sp.n2 = @(t) [ cos(sp.alpha + pi/2); sin(sp.alpha + pi/2) ]; % Second sliding sleeve normal vector

sp.p1_t = @(t) [0; 0];  % First and second time derivatives
sp.p1_tt = @(t) [0; 0];
sp.b1_t = @(t) [0; 0];
sp.b1_tt = @(t) [0; 0];
sp.n1_t = @(t) [0; 0];
sp.n1_tt = @(t) [0; 0];

sp.p2_t = @(t) [0; 0];
sp.p2_tt = @(t) [0; 0];
sp.b2_t = @(t) [0; 0];
sp.b2_tt = @(t) [0; 0];
sp.n2_t = @(t) [0; 0];
sp.n2_tt = @(t) [0; 0];


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

s1 = sp.s10;
v1 = sp.v10;
s2 = sp.s20;
v2 = sp.v20;

a1 = 0;
a2 = 0;

p = [s1; s2];
pp = [v1;v2];
ppp= [a1;a2];
p_prev = p;
pp_prev = pp;
ppp_prev = ppp;

fig = figure;
hold on
plot([U(sp.map.X(:,1)); U(sp.map.X(end,5))], [U(sp.map.X(:,2)); U(sp.map.X(end,6))], '-d' )
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
        [UM, F, J11] = newmark_matrices_two_sleeves( U, U_prev, Up_prev, Upp_prev, p, p_prev, pp_prev, ppp_prev, t, t_prev, sp );
        Jp1 = jacobian_numerical( @(U_var) interface_condition_two_sleeves(U_var, U_prev, Up_prev, Upp_prev, p, p_prev, pp_prev, ppp_prev, t, t_prev, sp), U, 1e-4);
        Jp  = jacobian_numerical( @(p_var) nl_system_extended(U, U_prev, Up_prev, Upp_prev, p_var, p_prev, pp_prev, ppp_prev, t, t_prev, sp), p, 1e-4);
        J1p = Jp(1:(5*(sp.Nel+1)+6),:);
        Jpp = Jp((5*(sp.Nel+1)+6)+1:end,:);

        FF = nl_system_extended( U, U_prev, Up_prev, Upp_prev, p, p_prev, pp_prev, ppp_prev, t, t_prev, sp );
        
        Jxx = J11(1:(4*(sp.Nel+1)), 1:(4*(sp.Nel+1)));
        Jxc = J11(1:(4*(sp.Nel+1)), 4*(sp.Nel+1) + [1:(sp.Nel+1+6)]);
        Jcx = J11(4*(sp.Nel+1) + [1:(sp.Nel+1+6)], 1:4*(sp.Nel+1));
        Jcc = J11(4*(sp.Nel+1) + [1:sp.Nel+1+6], 4*(sp.Nel+1) + [1:sp.Nel+1+6]);
        
        [Pe,Re,Ce] = equilibrate(Jxx);
        B = Re*Pe*Jxx*Ce;
        Jxxinv = Ce*( B\(Re*Pe) );

        J11sc  = Jcc - Jcx*Jxxinv*Jxc;
        J11sci = inv(J11sc);
        J11inv = [ Jxxinv + Jxxinv*Jxc*J11sci*Jcx*Jxxinv,  -Jxxinv*Jxc*J11sci;  -J11sci*Jcx*Jxxinv,  J11sci ];
        if cnt < 3
            J11inv = inv(J11);
        end
        
        Jsc = Jpp - Jp1*J11inv*J1p;
        Jsci = inv(Jsc);
        Jinv11 = J11inv + J11inv*J1p*Jsci*Jp1*J11inv;
        Jinv12 = -J11inv*J1p*Jsci;
        Jinv21 = -Jsci*Jp1*J11inv;
        Jinv22 = Jsci;
        
        U = U - ( Jinv11*FF(1:end-length(p)) + Jinv12*FF(end-length(p)+1:end) );
        p = p - ( Jinv21*FF(1:end-length(p)) + Jinv22*FF(end-length(p)+1:end) );
        
        residual = norm(UM * U - F);
        step_size = norm(U-U_old);
        step_dir = (U-U_old); step_dir = step_dir/max(abs(step_dir));
        disp([num2str(kk),': res=', num2str(residual),' ,  step=', num2str(step_size)]);
        if residual < 1e-9*sp.map.Ndof && step_size < 1e-9*sp.map.Ndof
            break;
        end
    end
    
    Upp = 1/(dt*dt*sp.newmark_b1)*(U-U_prev) - 1/(dt*sp.newmark_b1)*Up_prev - (1-2*sp.newmark_b1)/(2*sp.newmark_b1)*Upp_prev;
    Up = Up_prev + (1-sp.newmark_b2)*dt*Upp_prev + sp.newmark_b2*dt*Upp;

    ppp = 1/(dt*dt*sp.newmark_b1)*(p-p_prev) - 1/(dt*sp.newmark_b1)*pp_prev - (1-2*sp.newmark_b1)/(2*sp.newmark_b1)*ppp_prev; 
    pp  = pp_prev + (1-sp.newmark_b2)*dt*ppp_prev + sp.newmark_b2*dt*ppp;
    
    time(cnt) = t;
    DOFS(:, cnt) = [U; p];
    DOFS_p(:,cnt) = [Up; pp];
    DOFS_pp(:,cnt) = [Upp; ppp];

    if mod(cnt,20)==0
        clf(fig);
        subplot(2,1, 1)
        hold on
        plot([U(sp.map.X(:,1)); U(sp.map.X(end,5))], [U(sp.map.X(:,2)); U(sp.map.X(end,6))], '-o' )
        axis equal

        subplot(2,1,2)
        hold on
        plot(time, DOFS(end-1,:))
        plot(time, DOFS(end,:))
        yline(sp.L, '--k')
        yline(0, '--k')

        drawnow
    end

    if p(1)<0 || p(2)>sp.L
        break;
    end

    U_prev = U;
    Up_prev = Up;
    Upp_prev = Upp;

    p_prev = p;
    pp_prev = pp;
    ppp_prev = ppp;
end


function Fval = nl_system_extended( U, U_prev, Up_prev, Upp_prev, p, p_prev, pp_prev, ppp_prev, t, t_prev, sp )
    [UM, F] = newmark_matrices_two_sleeves( U, U_prev, Up_prev, Upp_prev, p, p_prev, pp_prev, ppp_prev, t, t_prev, sp );
    IC = interface_condition_two_sleeves(U, U_prev, Up_prev, Upp_prev, p, p_prev, pp_prev, ppp_prev, t, t_prev, sp);
    Fval = [UM*U - F; IC];
end