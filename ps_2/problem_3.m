% Define Q as an anonymous function
Q = @(x) (0.*(x>=0 & x<1) + ...
          abs(sin(6*(x - 1))).*(x>=1 & x<2) + ...
          (sqrt(x - 2) - sin(6)).*(x>=2 & x<3) + ...
          0.*(x>=3 & x<=4));

% Parameters for the problem
A = 9; % Z(0) = A
B = 8; % Z'(4) = B
alpha = sqrt(2);
N_fd = 100; % Number of points for the finite difference method
a = 0; % Start of interval
b = 4; % End of interval
dx_fd = (b - a) / (N_fd-1);
xx_fd = linspace(a, b, N_fd);
N_series = 40; % 

Q_values_fd = Q(xx_fd);

% Coefficients for the tridiagonal matrix
main_diag = (-2 + dx_fd^2 * alpha^2) * ones(N_fd, 1);
off_diag = ones(N_fd - 1, 1);
rhs = -dx_fd^2 * Q_values_fd;

% Solve the tridiagonal system using the tridiag_dB function
yy_fd = tridiag_dB([0; off_diag], main_diag, [off_diag; 0], rhs, A, B);

% Plot the finite difference solution
plot(xx_fd, yy_fd, 'b', 'LineWidth', 2); hold on;

% Call the Qfunc to plot the series solution, passing Q as an argument
Qfunc(N_fd, Q);

% Add the Wronskian method plot
%yy_wronskian = wronskianMethod(xx_fd);
%plot(xx_fd, yy_wronskian, 'g', 'LineWidth', 2);

% Adjust the figure properties
title('Solutions Comparison');
xlabel('x');
ylabel('Z(x)');
legend('Finite Difference', 'Series Expansion');
grid on;
hold off; % Stop adding to the current figure


% Function for the tridiagonal matrix solution
function yy = tridiag_dB(a,b,c,d,A,B)
    N = length(b);
    c(1) = c(1) / b(1);
    d(1) = (d(1) - A) / b(1);

    for i = 2:N
        temp = b(i) - a(i) * c(i-1);
        c(i) = c(i) / temp;
        d(i) = (d(i) - a(i) * d(i-1)) / temp;
    end

    yy = zeros(N, 1);
    yy(N) = (d(N) - B) / b(N);

    for i = N-1:-1:1
        yy(i) = d(i) - c(i) * yy(i+1);
    end
end

% Function for the series expansion solution
function Qfunc(N, Q)
    dx = 1/100;
    x = 0:dx:4;
    lam = ((2*(1:N)-1)*pi)/8;
    w = ones(size(x));
    q = Q(x);
    
    Qfunc = zeros(size(x));
    for n = 1:N
        yn = sin(lam(n)*x);
        Nn = simpson(yn.^2 .* w, dx);
        Q_hat = simpson(yn.*w, dx);
        Z_hat = Q_hat./(2-lam(n).^2);
        Qfunc = Qfunc + Z_hat * yn / Nn;
    end
    
    plot(x, Qfunc, '--r', 'LineWidth', 2);
end

% Function for Simpson's rule (numerical integration)
function foo = simpson(Q, dx)
    foo = dx/3 * (Q(1) + Q(end) + 4 * sum(Q(2:2:end-1)) + 2 * sum(Q(3:2:end-2)));
end

% Function for the Wronskian method solution
%function yy = wronskianMethod(x)
    %yp = @(x) ; %
    %C1 = 
    %C2 =
  
    %y_general = @(x) C1 .* x.^(-1) + C2 .* x + yp(x);
    
    % Evaluate the general solution
    %yy = y_general(x);
%end
