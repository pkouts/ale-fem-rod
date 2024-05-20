function Fval = interface_condition(U, U_prev, Up_prev, Upp_prev, p, p_prev, pp_prev, ppp_prev, t, t_prev, sp)
    map = sp.map;
    
    dt = t-t_prev;

    beta1 = sp.newmark_b1;
    beta2 = sp.newmark_b2;

    
    Upp = 1/(dt*dt*sp.newmark_b1)*(U-U_prev) - 1/(dt*sp.newmark_b1)*Up_prev - (1-2*sp.newmark_b1)/(2*sp.newmark_b1)*Upp_prev;
    Up = Up_prev + (1-beta2)*dt*Upp_prev + beta2*dt*Upp;
    
    ppp = 1/(dt*dt*sp.newmark_b1)*(p-p_prev) - 1/(dt*sp.newmark_b1)*pp_prev - (1-2*sp.newmark_b1)/(2*sp.newmark_b1)*ppp_prev; 
    pp  = pp_prev + (1-beta2)*dt*ppp_prev + beta2*dt*ppp; 
    % a1 = 1/(dt*dt*sp.newmark_b1)*(s1-s1_prev) - 1/(dt*sp.newmark_b1)*v1_prev - (1-2*sp.newmark_b1)/(2*sp.newmark_b1)*a1_prev;
    % v1 = v1_prev + (1-beta2)*dt*a1_prev + beta2*dt*a1;
    
    s1 = p(1); s1_prev = p_prev(1); v1_prev = pp_prev(1); a1_prev = ppp_prev(1); v1 = pp(1); a1 = ppp(1);
    s2 = p(2); s2_prev = p_prev(2); v2_prev = pp_prev(2); a2_prev = ppp_prev(2); v2 = pp(2); a2 = ppp(2);


    B = sp.B;
    % m = sp.m;
    gamma = sp.gamma;
    L = sp.L;
    g = sp.g;
    % ze = sp.ze;
    mu = sp.mu;

    s_node = mesh_to_material(sp.nodes, s1, s2, sp);

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

    R = U(map.clamp1(1:2));
    M = U(map.clamp1(3));
    p    = sp.p1(t);
    p_t  = sp.p1_t(t);
    p_tt = sp.p1_tt(t);
    b    = sp.b1(t);
    b_t  = sp.b1_t(t);
    b_tt = sp.b1_tt(t);
    n    = sp.n1(t);
    n_t  = sp.n1_t(t);
    n_tt = sp.n1_tt(t);

    epsilon = 1e-6;
    Ffriction = mu * sqrt(R'*n*n'*R) * v1/sqrt(v1*v1+epsilon);

    TC_s1 = (b*v1+b_t*s1-p_t)' * (b*v1+b_t*s1-p_t)*(gamma*0.5);
    VC_s1 = g'*( p - b*s1 ) * (-gamma);
    TC_s1_t_t = ( (b_t*s1+b*v1)' * (b*2*v1 + b_t*s1 - p_t*2) + b'*( (b_t*v1+b*a1)*2 + b_tt*s1 + b_t*v1 - p_tt*2 )*s1 ) * (gamma*0.5);

    

    Fval(1,:) = ( (Vs1 - Xds1*v1)'*(Vs1 - Xds1*v1)*(gamma*0.5) - Xdds1'*Xdds1*(B*0.5) + g'*Xs1*gamma )*(-1) + TC_s1 - TC_s1_t_t - VC_s1 - R'*Xds1 - M*n'*Xdds1 - Ffriction;

    % Fval(1,:) = s1 - ( sp.s10  );






    i = sp.Nel;

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


    Hs2   = interp_mat_C1_2D( 1, ds );
    Hds2  = interp_mat_C1_2D_der1( 1, ds ) / ds;
    Hdds2 = interp_mat_C1_2D_der2( 1, ds ) / (ds*ds);
    Xs2   = Hs2   * Xi;
    Xds2  = Hds2  * Xi;
    Xdds2 = Hdds2 * Xi;
    Vs2   = Hs2 * Vi;

    R2 = U(map.clamp2(1:2));
    M2 = U(map.clamp2(3));
    p2    = sp.p2(t);
    p2_t  = sp.p2_t(t);
    p2_tt = sp.p2_tt(t);
    b2    = sp.b2(t);
    b2_t  = sp.b2_t(t);
    b2_tt = sp.b2_tt(t);
    n2    = sp.n2(t);
    n2_t  = sp.n2_t(t);
    n2_tt = sp.n2_tt(t);

    epsilon = 1e-6;
    Ffriction = mu * sqrt(R2'*n2*n2'*R2) * v2/sqrt(v2*v2+epsilon);

    % VC = - ( gamma*(sp.L-s2)*g' ) * ( p2 + b2*(sp.L-s2)/2 );
    % TC = 0.5*gamma*(sp.L-s2)*(v2^2);

    % VC_s2 = - ( gamma*(-1)*g' ) * ( p2 + b2*(sp.L-s2)/2 ) - ( gamma*(sp.L-s2)*g' ) * ( b2*(-1)/2 );
    % TC_s2 = 0.5*gamma*(-1)*(v2^2);
    % TC_s2_t = gamma*(sp.L-s2)*(v2);
    % TC_s2_t_t = gamma*(-v2)*(v2) + gamma*(sp.L-s2)*(a2);


    VC_s2 = g'*( p2 + b2*(sp.L-s2) ) * (gamma);
    TC_s2 = ( b2_t*(sp.L-s2) + p2_t - b2*v2 )' * ( b2_t*(sp.L-s2) + p2_t - b2*v2 ) * (-gamma*0.5);
    
    A_s2_t = b2'*b2*(2*v2) - p2_t'*b2*2 + b2'*b2_t*(2*s2);
    B_s2_t = b2'*b2_t*(-2);

    A_s2_t_t = b2'*b2_t*(4*v2) + b2'*b2*(2*a2) - p2_tt'*b2*2 - p2_t'*b2_t*2 + b2'*b2_t*(2*v2) + b2_t'*b2_t*(2*s2) + b2'*b2_tt*(2*s2);
    B_s2_t_t = b2_t'*b2_t*(-2) + b2'*b2_tt*(-2);

    TC_s2_t_t = (  A_s2_t_t*(sp.L-s2) + A_s2_t*(-v2) + B_s2_t_t*(sp.L*sp.L-s2*s2)*0.5 + B_s2_t*(-s2*v2)  )*(0.5*gamma);

    % Fval(2,:) = s2 - sp.s20;
    Fval(2,:) = ( (Vs2 - Xds2*v2)'*(Vs2 - Xds2*v2)*(gamma*0.5) - Xdds2'*Xdds2*(B*0.5) + g'*Xs2*gamma )*(1) + TC_s2 - TC_s2_t_t - VC_s2 - R2'*Xds2 - n2'*Xdds2*M2 - Ffriction;


end