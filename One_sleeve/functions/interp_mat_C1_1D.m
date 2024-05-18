function mat = interp_mat_C1_1D( x, dx )
    mat = zeros(1,4);

    mat(1,1) = 1 + 0*x - 3*x*x + 2*x*x*x;
    mat(1,2) = (0 + 1*x - 2*x*x + 1*x*x*x)*dx;
    mat(1,3) = 0 + 0*x + 3*x*x - 2*x*x*x;
    mat(1,4) = (0 + 0*x - 1*x*x + 1*x*x*x)*dx;
end