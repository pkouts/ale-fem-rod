function quad = GaussQuadrature7()

    abscissas = zeros(7,1);
    weights = zeros(7,1);

    abscissas(1)  = 0.025446043828621;
    abscissas(2)  = 0.129234407200303;
    abscissas(3)  = 0.297077424311301;
    abscissas(4)  = 0.500000000000000;
    abscissas(5)  = 0.702922575688699;
    abscissas(6)  = 0.870765592799697;
    abscissas(7)  = 0.974553956171379;

    weights(1)  = 0.064742483084435;
    weights(2)  = 0.139852695744638;
    weights(3)  = 0.190915025252559;
    weights(4)  = 0.208979591836735;
    weights(5)  = 0.190915025252559;
    weights(6)  = 0.139852695744638;
    weights(7)  = 0.064742483084435;

    quad.abscissas = abscissas;
    quad.weights = weights;

end