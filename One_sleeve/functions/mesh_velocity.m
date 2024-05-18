function v = mesh_velocity( s, s1, v1, sp )
    v = v1*( 1-(s-s1)/(sp.L-s1) );
end