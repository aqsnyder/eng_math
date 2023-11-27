% Constants C1 and C2 from the solved initial conditions
C1 = 14.17;
C2 = -15.27;

% Define the particular solution yp and the general solution
yp = @(x) x.^3 .* (6*x + 5) / 10;
y_general = @(x) C1 * x.^(-1) + C2 * x + yp(x);

% Define the range for x
x = linspace(1, 5, 100);

% Evaluate the general solution over the range of x
y = y_general(x);

% Calculate the area under the curve using Simpson's rule
area_simpsons = simpsons_rule(y_general, 1, 5, 100);
disp(['Area under the curve using Simpson''s rule: ', num2str(area_simpsons)]);

% Plotting
plot(x, y, 'DisplayName', 'General Solution');
hold on;
title('General Solution of the Differential Equation');
xlabel('x');
ylabel('y');
legend show;
grid on;

% Simpson's rule function
function area = simpsons_rule(f, a, b, n)
    if mod(n, 2) == 1
        error('n must be an even integer.');
    end
    h = (b - a) / n;
    x = linspace(a, b, n+1);
    y = f(x);
    S = h/3 * sum(y(1:2:end-2) + 4*y(2:2:end-1) + y(3:2:end));
    area = S;
end
