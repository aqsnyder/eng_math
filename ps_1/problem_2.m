% Number of terms for power series
N = 50;

% Coefficient generation for power series
a = zeros(1, N+1);
a(1) = 1;
a(2) = 0;

for p = 0:N-1
    a(p+3) = -9*a(p+1)/((p+2)*(p+1));
end

% Arbitrary constants for the exact solution
C1 = 1;
C2 = 0;

% Values for x
x_vals = linspace(1, 5, 400);
y_approx_vals = arrayfun(@(x) y_approx(x, a, N), x_vals);
y_exact_vals = arrayfun(@(x) y_exact(x, C1, C2), x_vals);

% Plot
figure;
plot(x_vals, y_approx_vals, 'LineWidth', 2, 'Color', 'blue', 'DisplayName', 'Power Series Approximation');
hold on;
plot(x_vals, y_exact_vals, 'LineStyle', '--', 'LineWidth', 2, 'Color', 'red', 'DisplayName', 'Exact Solution');
title('Comparison of Power Series Approximation and Exact Solution');
xlabel('x');
ylabel('y(x)');
legend('show', 'Location', 'southeast');
grid on;

% Power series approximation function
function y = y_approx(x, a, N)
    y = 0;
    for p = 0:N
        y = y + a(p+1) * (x^p);
    end
end

% Exact solution function
function y = y_exact(x, C1, C2)
    y = C1 * cos(3*x) + C2 * sin(3*x);
end
