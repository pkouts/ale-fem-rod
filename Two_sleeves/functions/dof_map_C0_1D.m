function map = dof_map_C0_1D( Nel )

    map = zeros(Nel,2);

    for i=0:Nel-1
        map(i+1,1) = i;
        map(i+1,2) = i+1;
    end

end