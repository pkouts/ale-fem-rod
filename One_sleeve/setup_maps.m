function map = setup_maps(Nel)

    offset = 1;
    map.X = dof_map_C1_2D(Nel) + offset;
    offset = offset + (Nel+1)*4;
    map.N = dof_map_C0_1D(Nel) + offset;
    offset = offset + (Nel+1);
    map.clamp = dof_map_clamp_2D() + offset;
    offset = offset + 3;
    % map.s1 = offset;
    % offset = offset + 1;
    map.Ndof = offset-1;

end