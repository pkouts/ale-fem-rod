function s = mesh_to_material( xi, s1, s2, sp)
    s = xi*(s2 - s1) + s1;
end