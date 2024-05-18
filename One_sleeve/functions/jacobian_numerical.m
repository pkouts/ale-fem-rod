function J = jacobian_numerical( f, U, epsilon )
    f0 = f(U);
    J = zeros(length(f0), length(U));

    % epsilon = 1e-3;
    for ii=1:length(U)
        Up = U; Up(ii) = U(ii) + epsilon;
        Um = U; Um(ii) = U(ii) - epsilon;

        J(:,ii) = (f(Up) - f(Um)) / (2*epsilon);
    end
end