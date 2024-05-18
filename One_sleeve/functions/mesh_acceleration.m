function a = mesh_acceleration( s, s1, a1, sp )
    a = a1*( 1-(s-s1)/(sp.L-s1) );
end