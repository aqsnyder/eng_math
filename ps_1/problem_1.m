function main
    t = linspace(0, 10, 5000);
    dt = t(2) - t(1);
    
    v = zeros(size(t));
    v(1) = 0;

    for i = 2:length(t)
        v(i) = rk4(v(i-1), t(i-1), dt);
    end

    % Handle for v=0.95
    idx_95 = find(v >= 0.95, 1, 'first');
    if ~isempty(idx_95)
        time_95 = t(idx_95);
        fprintf('Velocity reaches 95%% at %.2f non-dimensional units.\n', time_95);
    else
        idx_95 = [];
        time_95 = [];
    end

    % Handle for v=1 (terminal velocity)
    idx_1 = find(v >= 0.999999, 1, 'first');
    if ~isempty(idx_1)
        time_1 = t(idx_1);
        fprintf('Velocity reaches terminal velocity at %.2f non-dimensional units.\n', time_1);
    else
        idx_1 = [];
        time_1 = [];
    end

    figure('Position', [100, 100, 800, 480]);
    plot(t, v, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Dimensionless Velocity $\hat{V}$');
    xlabel('Time (non-dimensional units)', 'Interpreter', 'latex');
    ylabel('Dimensionless Velocity $\hat{V}$', 'Interpreter', 'latex');
    title('Dimensionless Velocity vs Time', 'Interpreter', 'latex');
    grid on;
    hold on;

    if ~isempty(time_95)
        yline(0.95, 'r--', 'DisplayName', '95% of Terminal Velocity');
        xline(time_95, 'g--', 'DisplayName', sprintf('Time (95%%) = %.2f non-dimensional units', time_95));
    end

    if ~isempty(time_1)
        yline(1.0, 'm--', 'DisplayName', 'Terminal Velocity');
        xline(time_1, 'Color', [1 0.6 0], 'LineStyle', '--', 'DisplayName', sprintf('Time (Terminal) = %.2f non-dimensional units', time_1));
    end

    legend('show', 'Location', 'best', 'Interpreter', 'latex');
    hold off;
end

function dvdt = model(v)
    dvdt = 1 - v.^2;
end

function v_new = rk4(v, t, dt)
    k1 = dt * model(v);
    k2 = dt * model(v + 0.5 * k1);
    k3 = dt * model(v + 0.5 * k2);
    k4 = dt * model(v + k3);
    
    v_new = v + (k1 + 2*k2 + 2*k3 + k4) / 6.0;
end