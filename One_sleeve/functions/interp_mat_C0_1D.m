function mat = interp_mat_C0_1D( x )
    mat = zeros(1,2);

    mat(1,1) = 1-x;
    mat(1,2) = x;
end