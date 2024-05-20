function mat = interp_mat_C1_1D_der2( x, dx )
    mat = zeros(1,4);

    mat(1,1) = - 6 + 12*x;
    mat(1,2) = (- 4 + 6*x)*dx;
    mat(1,3) =   6 - 12*x;
    mat(1,4) = (- 2 + 6*x)*dx;
end