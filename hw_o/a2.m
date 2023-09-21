function cooling_simulation()

% Parameters
k = 0.1;  % Cooling constant
Ta = 25;  % Ambient temperature (degrees Celsius)

% Initial conditions
T0 = 100;  % Initial temperature (degrees Celsius)
t0 = 0;    % Initial time
tf = 10;   % Final time

% Time step and number of steps
dt = 0.1;
num_steps = (tf - t0) / dt;

% Arrays to store results
time_euler = zeros(1, num_steps + 1);
temp_euler = zeros(1, num_steps + 1);
time_predictor_corrector = zeros(1, num_steps + 1);
temp_predictor_corrector = zeros(1, num_steps + 1);

% Euler's method
time_euler(1) = t0;
temp_euler(1) = T0;
for i = 1:num_steps
    time_euler(i + 1) = time_euler(i) + dt;
    temp_euler(i + 1) = temp_euler(i) - k * (temp_euler(i) - Ta) * dt;
end

% Predictor-Corrector (Improved Euler) method
time_predictor_corrector(1) = t0;
temp_predictor_corrector(1) = T0;
for i = 1:num_steps
    time_predictor_corrector(i + 1) = time_predictor_corrector(i) + dt;
    % Predictor step
    predictor_temp = temp_predictor_corrector(i) - k * (temp_predictor_corrector(i) - Ta) * dt;
    % Corrector step
    temp_predictor_corrector(i + 1) = temp_predictor_corrector(i) - 0.5 * k * ((temp_predictor_corrector(i) - Ta) + (predictor_temp - Ta)) * dt;
end

% Plot results
figure;
plot(time_euler, temp_euler, 'DisplayName', "Euler's Method");
hold on;
plot(time_predictor_corrector, temp_predictor_corrector, 'DisplayName', "Predictor-Corrector (Improved Euler)");
xlabel('Time');
ylabel('Temperature (°C)');
title('Cooling of an Object');
legend;
grid on;
hold off;

end
