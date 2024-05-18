function map = dof_map_clamp_2D()

    map = zeros(1,3);

    map(1,1) = 0;
    map(1,2) = 1;
    map(1,3) = 2;

end