function [UM, F, J] = newmark_matrices( U, U_prev, Up_prev, Upp_prev, p, p_prev, pp_prev, ppp_prev, t, t_prev, sp )
    UM  = zeros(sp.map.Ndof);
    F   = zeros(sp.map.Ndof, 1);

    J = zeros(sp.map.Ndof);

    % Jxx2 = zeros(4*(sp.Nel+1));

    map = sp.map;
    
    dt = t-t_prev;

    beta1 = sp.newmark_b1;
    beta2 = sp.newmark_b2;

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
    % mu = sp.mu;
    c_distr = sp.c_distr;

    s_node = mesh_to_material(sp.nodes, s1, s2, sp);

    quadrature = GaussQuadrature7();

    for i=1:sp.Nel
        dxi = sp.nodes(i+1) - sp.nodes(i);
        ds  = s_node(i+1) - s_node(i);

        Xi = U( map.X(i,:) );
        Ni = U( map.N(i,:) );

        Xi_prev = U_prev( map.X(i,:) );
        Vi_prev = Up_prev( map.X(i,:) );
        Ai_prev = Upp_prev( map.X(i,:) );

        for k=1:length(quadrature.abscissas)
            x = quadrature.abscissas(k);
            w = quadrature.weights(k);

            H   = interp_mat_C1_2D( x, ds );
            Hd  = interp_mat_C1_2D_der1( x, ds ) / ds;
            Hdd = interp_mat_C1_2D_der2( x, ds ) / (ds*ds);
            P   = interp_mat_C0_1D( x );

            s = mesh_to_material( sp.nodes(i)+x*dxi, s1, s2, sp );
            v = mesh_velocity( s, s1, s2, v1, v2, sp );
            a = mesh_acceleration( s, s1, s2, a1, a2, sp );

            X   = H   * Xi;
            Xd  = Hd  * Xi;
            Xdd = Hdd * Xi;
            N   = P   * Ni;


            Mxx = ( H'*H*gamma ) * (ds*w);
            Kxx = ( Hd'*Hd * ( -gamma*v*v ) + Hdd'*Hdd*B - H'*Hd*a*gamma ) * (ds*w);
            Kxn = ( 2*Hd'*Xd*P ) * (ds*w);
            Knx = ( P'*Xd'*Hd ) * (ds*w);
            Cxx = ( H'*H*(-v1)/(L-s1)*gamma + Hd'*H*gamma*v - v*gamma*H'*Hd +  c_distr * H'*(eye(2) - Xd*Xd')*H  ) * (ds*w);

            Fx = ( H'*g*gamma ) * (ds*w);
            Fn = ( P' ) * (ds*w);

            Uxx = Kxx + Mxx*(1/(dt*dt*beta1)) + Cxx*(beta2/(dt*beta1));
            Uxn = Kxn;
            Unx = Knx;

            bx = Fx + Mxx*( (1/(dt*dt*beta1))*Xi_prev + (1/(dt*beta1))*Vi_prev + (1-2*beta1)/(2*beta1)*Ai_prev ) - Cxx*(Vi_prev + dt*(1-beta2-beta2*(1-2*beta1)/(2*beta1))*Ai_prev - beta2/beta1*Vi_prev - beta2/(dt*beta1)*Xi_prev);
            bn = Fn;

            Jxx = Uxx + ( 2*Hd'*Hd*(P*Ni) ) * (ds*w);
            Jxn = Uxn;
            Jnx = Uxn';

            % Jxx2(map.X(i,:), map.X(i,:)) = Jxx2(map.X(i,:), map.X(i,:)) + ( 2*Hd'*Hd*(P*Ni) ) * (ds*w);

            UM(map.X(i,:), map.X(i,:)) = UM(map.X(i,:), map.X(i,:)) + Uxx;
            UM(map.X(i,:), map.N(i,:)) = UM(map.X(i,:), map.N(i,:)) + Uxn;
            UM(map.N(i,:), map.X(i,:)) = UM(map.N(i,:), map.X(i,:)) + Unx;
            
            F(map.X(i,:)) = F(map.X(i,:)) + bx;
            F(map.N(i,:)) = F(map.N(i,:)) + bn;

            J(map.X(i,:), map.X(i,:)) = J(map.X(i,:), map.X(i,:)) + Jxx;
            J(map.X(i,:), map.N(i,:)) = J(map.X(i,:), map.N(i,:)) + Jxn;
            J(map.N(i,:), map.X(i,:)) = J(map.N(i,:), map.X(i,:)) + Jnx;
            
        end

        if i==1
            Hs1   = interp_mat_C1_2D( 0, ds );
            Hds1  = interp_mat_C1_2D_der1( 0, ds ) / ds;
            Hdds1 = interp_mat_C1_2D_der2( 0, ds ) / (ds*ds);

            Xs1   = Hs1   * Xi;
            Xds1  = Hds1  * Xi;
            Xdds1 = Hdds1 * Xi;
            % Vs1   = Hs1 * Vi;

            R = U(map.clamp1(1:2));
            M = U(map.clamp1(3));
            % R_prev = U_prev(map.clamp(1:2));
            % M_prev = U_prev(map.clamp(3));
            % Rp_prev = Up_prev(map.clamp(1:2));
            % Mp_prev = Up_prev(map.clamp(3));
            % Rpp_prev = Upp_prev(map.clamp(1:2));
            % Mpp_prev = Upp_prev(map.clamp(3));

            p1 = sp.p1(t);
            b1 = sp.b1(t);
            n1 = sp.n1(t);

            Uxr = Hs1';
            Uxm = Hds1'*n1;
            Urx = Hs1;
            Umx = n1'*Hds1;

            br = p1;
            bm = 0;
            % bxs1 = 0;

            Jxr = Uxr;
            Jxm = Uxm;
            Jrx = Urx;
            Jmx = Umx;

            UM(map.X(i,:), map.clamp1(1:2)) = UM(map.X(i,:), map.clamp1(1:2)) + Uxr;
            UM(map.X(i,:), map.clamp1(3))   = UM(map.X(i,:), map.clamp1(3))   + Uxm;
            UM(map.clamp1(1:2), map.X(i,:)) = UM(map.clamp1(1:2), map.X(i,:)) + Urx;
            UM(map.clamp1(3), map.X(i,:))   = UM(map.clamp1(3), map.X(i,:))   + Umx;
            % F(map.X(i,:)) = F(map.X(i,:)) + bxs1;
            F(map.clamp1) = F(map.clamp1) + [ br; bm ];

            J(map.X(i,:), map.clamp1(1:2)) = J(map.X(i,:), map.clamp1(1:2)) + Jxr;
            J(map.X(i,:), map.clamp1(3))   = J(map.X(i,:), map.clamp1(3))   + Jxm;
            J(map.clamp1(1:2), map.X(i,:)) = J(map.clamp1(1:2), map.X(i,:)) + Jrx;
            J(map.clamp1(3), map.X(i,:))   = J(map.clamp1(3), map.X(i,:))   + Jmx;
        end

        if i==sp.Nel
            HL   = interp_mat_C1_2D( 1, ds );
            HdL  = interp_mat_C1_2D_der1( 1, ds ) / ds;

            R = U(map.clamp2(1:2));
            M = U(map.clamp2(3));

            p2 = sp.p2(t);
            b2 = sp.b2(t);
            n2 = sp.n2(t);

            Uxr = HL';
            Uxm = HdL'*n2;
            Urx = HL;
            Umx = n2'*HdL;
            
            br = p2;
            bm = 0;

            % XL = HL * Xi;
            % VL = HL * Vi;
            % AL = HL * Ai;

            % disip = 2*ze*sqrt(3*m*B/(L-s1)^3);
            
            % MxxL = HL'*HL*m;
            % CxxL = HL'*HL*(disip);
            % FxL  = HL'*(g*m) + HL'*(sp.F1(t));
            
            % UxxL = ( MxxL*(1/(dt*dt*beta1)) + CxxL*(beta2/(dt*beta1)) );
            % bxL = FxL + MxxL*( (1/(dt*dt*beta1))*Xi_prev + (1/(dt*beta1))*Vi_prev + (1-2*beta1)/(2*beta1)*Ai_prev ) - CxxL*( Vi_prev + dt*(1-beta2-beta2*(1-2*beta1)/(2*beta1))*Ai_prev - beta2/beta1*Vi_prev - beta2/(dt*beta1)*Xi_prev );


            % UxxL = MxxL*(1/(dt*dt*beta1));
            % bxL  = FxL + MxxL*( 1/(dt*dt*beta1)*Xi_prev + 1/(dt*beta1)*Vi_prev + (1-2*beta1)/(2*beta1)*Ai_prev );

            % JxxL = UxxL;

            Jxr = Uxr;
            Jxm = Uxm;
            Jrx = Urx;
            Jmx = Umx;

            UM(map.X(i,:), map.clamp2(1:2)) = UM(map.X(i,:), map.clamp2(1:2)) + Uxr;
            UM(map.X(i,:), map.clamp2(3))   = UM(map.X(i,:), map.clamp2(3))   + Uxm;
            UM(map.clamp2(1:2), map.X(i,:)) = UM(map.clamp2(1:2), map.X(i,:)) + Urx;
            UM(map.clamp2(3), map.X(i,:))   = UM(map.clamp2(3), map.X(i,:))   + Umx;
            % F(map.X(i,:)) = F(map.X(i,:)) + bxs1;
            F(map.clamp2) = F(map.clamp2) + [ br; bm ];

            J(map.X(i,:), map.clamp2(1:2)) = J(map.X(i,:), map.clamp2(1:2)) + Jxr;
            J(map.X(i,:), map.clamp2(3))   = J(map.X(i,:), map.clamp2(3))   + Jxm;
            J(map.clamp2(1:2), map.X(i,:)) = J(map.clamp2(1:2), map.X(i,:)) + Jrx;
            J(map.clamp2(3), map.X(i,:))   = J(map.clamp2(3), map.X(i,:))   + Jmx;
        end
    end

    % UM = UM * dt^(0) * ds^(2);
    % F  = F  * dt^(0) * ds^(2);
    % J  = J  * dt^(0) * ds^(2);
end
