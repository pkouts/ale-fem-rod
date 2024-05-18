function [U, Up, Upp] = initialize_dofs(t_init, sp)

    U  = zeros(sp.map.Ndof,1);
    Up = zeros(sp.map.Ndof,1);
    Upp = zeros(sp.map.Ndof,1);

    s1 = sp.L-sp.l0;
    v1 = -sp.l0_t;

    % U(sp.map.s1)  = s1;
    % Up(sp.map.s1) = v1;

    for i=1:sp.Nel
        s = mesh_to_material( sp.nodes(i), s1, sp);
        p = sp.p(t_init);
        p_t = sp.p_t(t_init);

        b = sp.b(t_init);
        b_t = sp.b_t(t_init);

        x = p + b*(s-s1);
        x_t = p_t - b*v1 + b_t*(s-s1);

        U(sp.map.X(i,1:4))  = [ x; b ];
        % Up(sp.map.X(i,1:4)) = [ x_t; b_t ];

        % U(sp.map.N(i,1)) = sp.m*sp.g'*b;

        if i==sp.Nel
            s = mesh_to_material( sp.nodes(i+1), s1, sp);
            x = p + b*(s-s1);
            x_t = p_t - b*v1 + b_t*(s-s1);
            U(sp.map.X(i,5:8))  = [ x; b ];
            % Up(sp.map.X(i,5:8)) = [ x_t; b_t ];
            % U(sp.map.N(i,2)) = sp.m*sp.g'*b;
        end
    end
end