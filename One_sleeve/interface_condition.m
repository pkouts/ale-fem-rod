function Fval = interface_condition(U, U_prev, Up_prev, Upp_prev, s1, s1_prev, v1_prev, a1_prev, t, t_prev, sp)
    map = sp.map;
    
    dt = t-t_prev;

    beta1 = sp.newmark_b1;
    beta2 = sp.newmark_b2;

    Upp = 1/(dt*dt*sp.newmark_b1)*(U-U_prev) - 1/(dt*sp.newmark_b1)*Up_prev - (1-2*sp.newmark_b1)/(2*sp.newmark_b1)*Upp_prev;
    % Up  = Up_prev + dt*( 1-sp.newmark_b2 - sp.newmark_b2*(1-2*sp.newmark_b1)/(2*sp.newmark_b1) )*Upp_prev - sp.newmark_b2/sp.newmark_b1*Up_prev - sp.newmark_b2/(dt*sp.newmark_b1)*U_prev + sp.newmark_b2/(dt*sp.newmark_b1)*U;
    Up = Up_prev + (1-beta2)*dt*Upp_prev + beta2*dt*Upp;

    a1 = 1/(dt*dt*sp.newmark_b1)*(s1-s1_prev) - 1/(dt*sp.newmark_b1)*v1_prev - (1-2*sp.newmark_b1)/(2*sp.newmark_b1)*a1_prev;
    % v1  = v1_prev + dt*( 1-sp.newmark_b2 - sp.newmark_b2*(1-2*sp.newmark_b1)/(2*sp.newmark_b1) )*a1_prev - sp.newmark_b2/sp.newmark_b1*v1_prev - sp.newmark_b2/(dt*sp.newmark_b1)*s1_prev + sp.newmark_b2/(dt*sp.newmark_b1)*s1;
    v1 = v1_prev + (1-beta2)*dt*a1_prev + beta2*dt*a1;

    B = sp.B;
    m = sp.m;
    gamma = sp.gamma;
    L = sp.L;
    g = sp.g;
    ze = sp.ze;
    mu = sp.mu;

    s_node = mesh_to_material(sp.nodes, s1, sp);

    i = 1;

    dxi = sp.nodes(i+1) - sp.nodes(i);
    ds  = s_node(i+1) - s_node(i);

    Xi = U( map.X(i,:) );
    Vi = Up( map.X(i,:) );
    Ai = Upp( map.X(i,:) );
    Ni = U( map.N(i,:) );

    Xi_prev = U_prev( map.X(i,:) );
    Vi_prev = Up_prev( map.X(i,:) );
    Ai_prev = Upp_prev( map.X(i,:) );
    Ni_prev = U_prev( map.N(i,:) );
    VNi_prev = Up_prev( map.N(i,:) );
    ANi_prev = Upp_prev( map.N(i,:) );


    Hs1   = interp_mat_C1_2D( 0, ds );
    Hds1  = interp_mat_C1_2D_der1( 0, ds ) / ds;
    Hdds1 = interp_mat_C1_2D_der2( 0, ds ) / (ds*ds);

    Xs1   = Hs1   * Xi;
    Xds1  = Hds1  * Xi;
    Xdds1 = Hdds1 * Xi;
    Vs1   = Hs1 * Vi;

    R = U(map.clamp(1:2));
    M = U(map.clamp(3));
    p    = sp.p(t);
    p_t  = sp.p_t(t);
    p_tt = sp.p_tt(t);
    b    = sp.b(t);
    b_t  = sp.b_t(t);
    b_tt = sp.b_tt(t);
    n    = sp.n(t);
    n_t  = sp.n_t(t);
    n_tt = sp.n_tt(t);

    epsilon = 1e-6;
    Ffriction = mu * sqrt(R'*n*n'*R) * v1/sqrt(v1*v1+epsilon);

    TC_s1 = (b*v1+b_t*s1-p_t)' * (b*v1+b_t*s1-p_t)*(gamma*0.5);
    VC_s1 = g'*( p - b*s1 ) * (-gamma);
    TC_s1_t_t = ( (b_t*s1+b*v1)' * (b*2*v1 + b_t*s1 - p_t*2) + b'*( (b_t*v1+b*a1)*2 + b_tt*s1 + b_t*v1 - p_tt*2 )*s1 ) * (gamma*0.5);

    % Fval = -0.5*B*Xdds1'*Xdds1 + R'*Xds1 + M*n'*Xdds1;

    Fval = ( (Vs1 - Xds1*v1)'*(Vs1 - Xds1*v1)*(gamma*0.5) - Xdds1'*Xdds1*(B*0.5) + g'*Xs1*gamma )*(-1) + TC_s1 - TC_s1_t_t - VC_s1 - R'*Xds1 - M*n'*Xdds1 - Ffriction;

    % Fval = gamma*b'*( (b*a1 - p_tt)*s1 + (b*v1 - p_t)*v1 + g*s1 ) - Xdds1'*Xdds1*(B*0.5) + R'*Xds1 + M*n'*Xdds1 + Ffriction;
    % Fval = TC_s1_t_t - TC_s1 + VC_s1 - Xdds1'*Xdds1*(B*0.5) + R'*Xds1 + M*n'*Xdds1 + Ffriction;

    % shilei_ell = readmatrix('heavy_oscillating_s1_shilei.csv');
    % s1_ref = interp1(shilei_ell(:,1), shilei_ell(:,2), t);
    % Fval = s1 - s1_ref;
    % Fval = s1-(sp.L-sp.l0);
    % fr = 5;
    % amp = 0.3;
    % Fval = s1 - (sp.L-sp.l0 + amp*(1-cos(t*2*pi*fr))).*(t<0.5/fr) - (sp.L-sp.l0 + amp*2).*(t>=0.5/fr);

    % Fval = Fval * dt^(0) * ds^(4);

end