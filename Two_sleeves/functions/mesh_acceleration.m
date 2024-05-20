function a = mesh_acceleration( s, s1, s2, a1, a2, sp )
    a = (s-s1)/(s2-s1) * (a2-a1) + a1;
end