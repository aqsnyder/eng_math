function problem_2(terms)
    weightFunction = @(x) x.^3;
    figure;
    for index = 1:length(terms)
        eigenvalues = computeEigenvalues(terms(index));
        numPoints = 4*length(eigenvalues); if numPoints<100, numPoints=100; end
        stepSize = 4/numPoints; xValues = 0:stepSize:4;
        targetFunction = computeQQ(xValues);
        rep = zeros(size(xValues));
        eigenFunction = zeros(size(xValues));

        for eigenvalue = eigenvalues
            eigenFunction(1) = eigenvalue/2;
            eigenFunction(2:end) = besselj(1, eigenvalue*xValues(2:end))./xValues(2:end);
            N_hat = simpsonIntegral(eigenFunction.^2.*weightFunction(xValues),stepSize);
            Q_hat = simpsonIntegral(eigenFunction.*targetFunction.*weightFunction(xValues),stepSize);
            rep = rep + Q_hat*eigenFunction/N_hat;
        end

        subplot(2, 2, index);
        plot(xValues, rep, ':', xValues, targetFunction, '-k', 'linewidth', 2);
        title(sprintf('Using %d terms of eigens', terms(index)));
    end
end

function integral = simpsonIntegral(values, dx)
    integral = dx/3*(values(1)+values(end)+4*sum(values(2:2:end-1))+2*sum(values(3:2:end-2)));
end

function eigenvalues = computeEigenvalues(terms)
    eigenvalues = zeros(1, terms); currentRootIndex = 0;
    functionValue = @(l) 0.5*(besselj(0, 4*l) - besselj(2, 4*l));
    x = 0; step = 1/10;
    while currentRootIndex < terms
        while functionValue(x) * functionValue(x+step) > 0, x = x + step; end
        currentRootIndex = currentRootIndex + 1;
        eigenvalues(currentRootIndex) = fzero(functionValue, x);
        x = x + 2*step;
    end
end

function QQvalues = computeQQ(x)
    % Nested function for QQ
    QQvalues = zeros(size(x));
    QQvalues(x >= 1 & x < 2) = abs(sin(6 * (x(x >= 1 & x < 2) - 1)));
    QQvalues(x >= 2 & x < 3) = sqrt(x(x >= 2 & x < 3) - 2) - sin(6);
    QQvalues(x >= 3) = 0;
end

% problem_2([5, 10, 20, 40]);
