function mat = interp_mat_C1_1D_der1( x, dx )
    mat = zeros(1,4);

    mat(1,1) = 0 - 6*x + 6*x*x;
    mat(1,2) = (1 - 4*x + 3*x*x)*dx;
    mat(1,3) = 0 + 6*x - 6*x*x;
    mat(1,4) = (0 - 2*x + 3*x*x)*dx;
end