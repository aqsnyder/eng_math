function main()

% Define the function to be integrated
f = @(x) 1 ./ (2 - sqrt(x));

% Integrate the function from 0 to 3.99 and from 4.01 to 5
integral_value_1 = integral(f, 0, 3.99999);
integral_value_2 = integral(f, 4.00001, 5);

% Sum the two integrals
total_integral = integral_value_1 + integral_value_2;

% Display the result
fprintf('Value of I from 0 to 5: %.5f\n', total_integral);

end
