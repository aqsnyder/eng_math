% Define the function v(t)
v = @(t) (exp(2*t) - 1) ./ (exp(2*t) + 1);

% Generate a range of t values
t = linspace(0, 5, 400);

% Calculate v for each t
v_values = v(t);

% Find the time t when velocity v is approximately 0.95
target_v = 0.95;
[~, idx] = min(abs(v_values - target_v));
t_at_target_v = t(idx);
fprintf('Non-dimensional time when velocity is approximately %.2f: %.2f\n', target_v, t_at_target_v);

% Create the plot
figure;
plot(t, v_values, 'DisplayName', 'v(t) = (e^{2t} - 1) / (e^{2t} + 1)');
title('Graph of v(t)');
xlabel('t');
ylabel('v(t)');
hold on;
grid on;

% Plotting lines
yline(0, 'k', 'LineWidth', 0.5);
xline(0, 'k', 'LineWidth', 0.5);
yline(target_v, 'r--', 'DisplayName', sprintf('v = %.2f', target_v));
xline(t_at_target_v, 'g--', 'DisplayName', sprintf('t at v = %.2f (%.2f)', target_v, t_at_target_v));

legend;
hold off;