function map = dof_map_C1_2D( Nel )

    map = zeros(Nel,8);

    for i=0:Nel-1
        map(i+1,1) = i;
        map(i+1,2) = (Nel+1)+i;
        map(i+1,3) = 2*(Nel+1)+i;
        map(i+1,4) = 3*(Nel+1)+i;
        map(i+1,5) = i+1;
        map(i+1,6) = (Nel+1)+i+1;
        map(i+1,7) = 2*(Nel+1)+i+1;
        map(i+1,8) = 3*(Nel+1)+i+1;
    end

end