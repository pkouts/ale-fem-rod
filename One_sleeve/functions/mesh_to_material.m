function s = mesh_to_material( xi, s1, sp)
    s = xi*(sp.L - s1) + s1;
end