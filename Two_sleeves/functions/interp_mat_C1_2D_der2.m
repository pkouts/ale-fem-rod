function mat = interp_mat_C1_2D_der2( x, dx )
    mat = zeros(2,8);
    h = interp_mat_C1_1D_der2(x,dx);

    mat(:,1:2) = eye(2)*h(1);
    mat(:,3:4) = eye(2)*h(2);
    mat(:,5:6) = eye(2)*h(3);
    mat(:,7:8) = eye(2)*h(4);

end