function v = mesh_velocity( s, s1, s2, v1, v2, sp )
    v = (s-s1)/(s2-s1) * (v2-v1) + v1;
end