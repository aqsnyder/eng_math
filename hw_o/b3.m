function runge_kutta_plot()

    TUID = [9, 1, 5, 1, 8, 7, 2, 8, 9];
    LETTER_MAP = {'I', 'H', 'G', 'F', 'E', 'D', 'C', 'B', 'A'};
    
    total_sum = sum(TUID);
    average = total_sum / length(TUID);
    
    fprintf('my TUID average: %.2f\n', average);
    
    % Create a map (in MATLAB, we use containers.Map) to map letters to integers
    letter_to_int_map = containers.Map(LETTER_MAP, TUID);
    
    % Calculate the averages of the letters of interest
    alpha_letters = {'A', 'B', 'C'};
    beta_letters = {'D', 'E', 'F'};
    gamma_letters = {'G', 'H', 'I'};
    
    alpha = mean(cellfun(@(x) letter_to_int_map(x), alpha_letters)) / 10;
    beta = mean(cellfun(@(x) letter_to_int_map(x), beta_letters)) / 10;
    gamma = mean(cellfun(@(x) letter_to_int_map(x), gamma_letters)) / 10;

    % Parameters
    T = 5;
    h = 0.01;

    % From TUID letter mapping
    A = 9;
    B = 8;
    C = 2;

    u0 = [A, B, C];
    
    [t_values, u_values] = runge_kutta_4th_order(h, T, u0, alpha, beta);
    
    figure;
    plot(t_values, u_values(1,:), 'DisplayName', 'y(t)');
    hold on;
    plot(t_values, u_values(2,:), 'DisplayName', "y'(t)");
    plot(t_values, u_values(3,:), 'DisplayName', "y''(t)");
    legend();
    xlabel('t');
    ylabel('Value');
    title('Integration using 4th order Runge-Kutta method');
    grid on;

    function [t_values, u_values] = runge_kutta_4th_order(h, T, u0, alpha, beta)
        % The system of ODEs
        f = @(t, u) [u(2); u(3); cos(3*t) - alpha*u(3) - beta*u(1)*u(3)];
        
        t_values = 0:h:T;
        u_values = zeros(3, length(t_values));
        u_values(:, 1) = u0;

        for i = 1:length(t_values)-1
            u = u_values(:, i);
            
            k1 = h * f(t_values(i), u);
            k2 = h * f(t_values(i) + 0.5*h, u + 0.5*k1);
            k3 = h * f(t_values(i) + 0.5*h, u + 0.5*k2);
            k4 = h * f(t_values(i) + h, u + k3);
            
            u_values(:, i+1) = u + (1/6) * (k1 + 2*k2 + 2*k3 + k4);
        end
    end
end
