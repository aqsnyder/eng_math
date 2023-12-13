% Parameters
alpha = 9;
beta = 8;
gamma = 5;
A = 0.5; % still need to solve for the actual coefficients in the Sturm-Liouville solution
% ran out of time for solving this problem 

% Finite Difference Method Setup
N = 60;
epsilon = 1e-3;
x = linspace(0, 2, N);
y = linspace(0, 1, N);
dx = x(2) - x(1);
dy = y(2) - y(1);
u = zeros(N, N);

% Apply boundary conditions
u(1, :) = beta * y;
u(end, :) = beta * y;

f = zeros(N, N);
for j = 1:N
    for i = round(N/2):N
        f(j, i) = 40 * sin(alpha * y(j)) * (1 - x(i));
    end
end

max_diff = inf;
iterations = 0;

while max_diff > epsilon
    max_diff = 0;
    for j = 2:(N-1)
        for i = 2:(N-1)
            u_old = u(j, i);
            if x(i) < 1
                u(j, i) = 0.25 * (u(j, i+1) + u(j, i-1) + u(j+1, i) + u(j-1, i));
            else
                u(j, i) = 0.25 * (u(j, i+1) + u(j, i-1) + u(j+1, i) + u(j-1, i) - dx^2 * f(j, i));
            end
            diff = abs(u(j, i) - u_old);
            max_diff = max(max_diff, diff);
        end
    end
    
    % Apply the Neumann boundary condition at y=1
    for i = 2:(N-1)
        u(end, i) = u(end-1, i) - dy * gamma * u(end-1, i);
    end
    
    iterations = iterations + 1;
end

% Sturm-Liouville solutions
k_root = 5.00045361; % Previously computed eigenvalue
sl_u_y_1 = A * cosh(k_root * y);
sl_u_x_2 = A * sin((pi/2) * x);

% Extracting finite difference solutions at the specified lines
x_index_1 = round(1.5 / dx);
y_index_2 = round(0.75 / dy);
fd_u_y_1 = u(:, x_index_1);
fd_u_x_2 = u(y_index_2, :);

% Plotting the comparison
figure;

% Plot for u(y) at x=3/2
subplot(1, 2, 1);
plot(y, fd_u_y_1, 'r--');
hold on;
plot(y, sl_u_y_1, 'k-');
xlabel('y');
ylabel('u');
title('u(y) at x=3/2 Comparison');
legend('Finite Difference Method', 'Sturm-Liouville Solution');
grid on;

% Plot for u(x) at y=3/4
subplot(1, 2, 2);
plot(x, fd_u_x_2, 'r--');
hold on;
plot(x, sl_u_x_2, 'k-');
xlabel('x');
ylabel('u');
title('u(x) at y=3/4 Comparison');
legend('Finite Difference Method', 'Sturm-Liouville Solution');
grid on;

% Surface plot of the solution u(x,y)
figure;
surf(x, y, u');
xlabel('X');
ylabel('Y');
zlabel('u');
title('Surface plot of the solution u(x,y)');
colorbar;
