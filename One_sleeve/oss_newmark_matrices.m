function [UM, F, J, Jxx2] = oss_newmark_matrices_alt( U, U_prev, Up_prev, Upp_prev, s1, s1_prev, v1_prev, a1_prev, t, t_prev, sp )
    UM  = zeros(sp.map.Ndof);
    F   = zeros(sp.map.Ndof, 1);

    J = zeros(sp.map.Ndof);

    Jxx2 = zeros(4*(sp.Nel+1));

    map = sp.map;
    
    dt = t-t_prev;

    beta1 = sp.newmark_b1;
    beta2 = sp.newmark_b2;

    % Upp = 1/(dt*dt*sp.newmark_b1)*(U-U_prev) - 1/(dt*sp.newmark_b1)*Up_prev - (1-2*sp.newmark_b1)/(2*sp.newmark_b1)*Upp_prev;
    % Up  = Up_prev + dt*( 1-sp.newmark_b2 - sp.newmark_b2*(1-2*sp.newmark_b1)/(2*sp.newmark_b1) )*Upp_prev - sp.newmark_b2/sp.newmark_b1*Up_prev - sp.newmark_b2/(dt*sp.newmark_b1)*U_prev + sp.newmark_b2/(dt*sp.newmark_b1)*U;
    % Up = Up_prev + (1-beta2)*dt*Upp_prev + beta2*dt*Upp;

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

            s = mesh_to_material( sp.nodes(i)+x*dxi, s1, sp );
            v = mesh_velocity( s, s1, v1, sp );
            a = mesh_acceleration( s, s1, a1, sp );

            X   = H   * Xi;
            Xd  = Hd  * Xi;
            Xdd = Hdd * Xi;
            N   = P   * Ni;

            % Mxx = ( H'*H*gamma ) * (ds*w);
            % Kxx = ( Hd'*Hd * ( -gamma*v*v ) + Hdd'*Hdd*B ) * (ds*w);
            % Kxn = ( 2*Hd'*Xd*P ) * (ds*w);
            % Knx = ( P'*Xd'*Hd ) * (ds*w);

            % Mxx = ( H'*H*gamma ) * (ds*w);
            % Kxx = ( Hd'*Hd * ( -gamma*v*v ) + Hdd'*Hdd*B ) * (ds*w);
            % Kxn = ( 2*Hd'*Xd*P ) * (ds*w);
            % Knx = ( P'*Xd'*Hd ) * (ds*w);
            % Cxx = ( zeros(size(H'*H)) ) * (ds*w);

            % Mxx = ( H'*H*gamma ) * (ds*w);
            % Kxx = ( Hd'*Hd * ( -gamma*v*v ) + Hdd'*Hdd*B - H'*Hd*a*gamma - H'*Hd*v*(-v1)/(L-s1)*gamma  ) * (ds*w);
            % Kxn = ( 2*Hd'*Xd*P ) * (ds*w);
            % Knx = ( P'*Xd'*Hd ) * (ds*w);
            % Cxx = ( H'*H*(-v1)/(L-s1)*gamma + Hd'*H*gamma*v ) * (ds*w);

            Mxx = ( H'*H*gamma ) * (ds*w);
            Kxx = ( Hd'*Hd * ( -gamma*v*v ) + Hdd'*Hdd*B - H'*Hd*a*gamma   ) * (ds*w);
            Kxn = ( 2*Hd'*Xd*P ) * (ds*w);
            Knx = ( P'*Xd'*Hd ) * (ds*w);
            Cxx = ( H'*H*(-v1)/(L-s1)*gamma + Hd'*H*gamma*v - v*gamma*H'*Hd ) * (ds*w);

            % Mxx = ( H'*H*gamma ) * (ds*w);
            % Kxx = ( Hd'*Hd * ( -gamma*v*v ) + Hdd'*Hdd*B - H'*Hd*a*gamma - H'*Hd*v*(-v1)/(L-s1)*gamma  ) * (ds*w);
            % Kxn = ( 2*Hd'*Xd*P ) * (ds*w);
            % Knx = ( P'*Xd'*Hd ) * (ds*w);
            % Cxx = ( -H'*H*(-v1)/(L-s1)*gamma ) * (ds*w);

            Fx = ( H'*g*gamma ) * (ds*w);
            Fn = ( P' ) * (ds*w);

            Uxx = Kxx + Mxx*(1/(dt*dt*beta1)) + Cxx*(beta2/(dt*beta1));
            Uxn = Kxn;
            Unx = Knx;

            bx = Fx + Mxx*( (1/(dt*dt*beta1))*Xi_prev + (1/(dt*beta1))*Vi_prev + (1-2*beta1)/(2*beta1)*Ai_prev ) - Cxx*(Vi_prev + dt*(1-beta2-beta2*(1-2*beta1)/(2*beta1))*Ai_prev - beta2/beta1*Vi_prev - beta2/(dt*beta1)*Xi_prev);
            bn = Fn;

            Jxx = Uxx + ( 2*Hd'*Hd*(P*Ni) ) * (ds*w);
            Jxn = Uxn;
            % Jnx = 2*Unx;
            Jnx = Uxn';

            Jxx2(map.X(i,:), map.X(i,:)) = Jxx2(map.X(i,:), map.X(i,:)) + ( 2*Hd'*Hd*(P*Ni) ) * (ds*w);

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

            R = U(map.clamp(1:2));
            M = U(map.clamp(3));
            % R_prev = U_prev(map.clamp(1:2));
            % M_prev = U_prev(map.clamp(3));
            % Rp_prev = Up_prev(map.clamp(1:2));
            % Mp_prev = Up_prev(map.clamp(3));
            % Rpp_prev = Upp_prev(map.clamp(1:2));
            % Mpp_prev = Upp_prev(map.clamp(3));

            p = sp.p(t);
            b = sp.b(t);
            n = sp.n(t);

            Uxr = Hs1';
            Uxm = Hds1'*n;
            Urx = Hs1;
            Umx = n'*Hds1;

            br = p;
            bm = 0;
            % bxs1 = 0;

            Jxr = Uxr;
            Jxm = Uxm;
            Jrx = Urx;
            Jmx = Umx;

            UM(map.X(i,:), map.clamp(1:2)) = UM(map.X(i,:), map.clamp(1:2)) + Uxr;
            UM(map.X(i,:), map.clamp(3))   = UM(map.X(i,:), map.clamp(3))   + Uxm;
            UM(map.clamp(1:2), map.X(i,:)) = UM(map.clamp(1:2), map.X(i,:)) + Urx;
            UM(map.clamp(3), map.X(i,:))   = UM(map.clamp(3), map.X(i,:))   + Umx;
            % F(map.X(i,:)) = F(map.X(i,:)) + bxs1;
            F(map.clamp) = F(map.clamp) + [ br; bm ];

            J(map.X(i,:), map.clamp(1:2)) = J(map.X(i,:), map.clamp(1:2)) + Jxr;
            J(map.X(i,:), map.clamp(3))   = J(map.X(i,:), map.clamp(3))   + Jxm;
            J(map.clamp(1:2), map.X(i,:)) = J(map.clamp(1:2), map.X(i,:)) + Jrx;
            J(map.clamp(3), map.X(i,:))   = J(map.clamp(3), map.X(i,:))   + Jmx;
        end

        if i==sp.Nel
            HL   = interp_mat_C1_2D( 1, ds );
            
            % XL = HL * Xi;
            % VL = HL * Vi;
            % AL = HL * Ai;

            disip = 2*ze*sqrt(3*m*B/(L-s1)^3);
            
            MxxL = HL'*HL*m;
            CxxL = HL'*HL*(disip);
            FxL  = HL'*(g*m) + HL'*(sp.F1(t));
            
            UxxL = ( MxxL*(1/(dt*dt*beta1)) + CxxL*(beta2/(dt*beta1)) );
            bxL = FxL + MxxL*( (1/(dt*dt*beta1))*Xi_prev + (1/(dt*beta1))*Vi_prev + (1-2*beta1)/(2*beta1)*Ai_prev ) - CxxL*( Vi_prev + dt*(1-beta2-beta2*(1-2*beta1)/(2*beta1))*Ai_prev - beta2/beta1*Vi_prev - beta2/(dt*beta1)*Xi_prev );


            % UxxL = MxxL*(1/(dt*dt*beta1));
            % bxL  = FxL + MxxL*( 1/(dt*dt*beta1)*Xi_prev + 1/(dt*beta1)*Vi_prev + (1-2*beta1)/(2*beta1)*Ai_prev );

            JxxL = UxxL;

            UM(map.X(i,:), map.X(i,:)) = UM(map.X(i,:), map.X(i,:)) + UxxL;
            F(map.X(i,:)) = F(map.X(i,:)) + bxL;

            J(map.X(i,:), map.X(i,:)) = J(map.X(i,:), map.X(i,:)) + JxxL;
        end
    end

    % UM = UM * dt^(0) * ds^(4);
    % F  = F  * dt^(0) * ds^(4);
    % J  = J  * dt^(0) * ds^(4);
end
