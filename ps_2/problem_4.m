% Set the number of grid points and the convergence criterion
N = 100;
alpha = 9;
beta = 8;
gamma = 5;
epsilon = 1e-3;

% Initialize the grid
u = zeros(N, N);
x = linspace(0, 2, N);
y = linspace(0, 1, N);
dx = x(2) - x(1);
dy = y(2) - y(1);

% Apply boundary conditions
u(1, :) = beta * y;  % u(x=0) = beta * y
u(N, :) = beta * y;  % u(x=2) = beta * y
% u(y=0) is already initialized to 0

% Precompute the constant f(x,y) for x >= 1
f = zeros(N, N);
for j = 1:N
    for i = round(N/2):N  % Start at x=1, which is halfway across the grid
        f(j, i) = 40 * sin(alpha * y(j)) * (1 - x(i));
    end
end

% Perform Gauss-Seidel iteration
max_diff = inf;
iterations = 0;

while max_diff > epsilon
    max_diff = 0;
    for j = 2:N-1
        for i = 2:N-1
            u_old = u(j, i);
            % Apply the Gauss-Seidel update
            if x(i) < 1
                u(j, i) = 0.25 * (u(j, i+1) + u(j, i-1) + u(j+1, i) + u(j-1, i));
            else
                u(j, i) = 0.25 * (u(j, i+1) + u(j, i-1) + u(j+1, i) + u(j-1, i) - dx^2 * f(j, i));
            end
            
            % Compute the change and track the maximum
            diff = abs(u(j, i) - u_old);
            max_diff = max(max_diff, diff);
        end
    end

    % Apply the Neumann boundary condition at y=1
    for i = 2:N-1
        u(N, i) = u(N-1, i) - dy * gamma * u(N-1, i);
    end
    
    iterations = iterations + 1;
end

% Create the meshgrid for plotting
[X, Y] = meshgrid(x, y);

% Plot the results
surf(X, Y, u', 'EdgeColor', 'none');
xlabel('X');
ylabel('Y');
zlabel('u');
title('Surface plot of the solution u(x,y)');
colorbar;

% Output the number of iterations and the maximum difference
disp(['Iterations: ', num2str(iterations)]);
disp(['Maximum change: ', num2str(max_diff)]);
