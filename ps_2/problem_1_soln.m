% List of Nterms to run
NtermsList = [5, 20, 40];

% Number of subplots
numSubplots = length(NtermsList);

% Create a figure for the subplots
figure;

% Execute problem_1 for each Nterms and place each in a subplot
for i = 1:numSubplots
    N = NtermsList(i);
    subplot(1, numSubplots, i); % Adjust the subplot grid as needed
    problem_1(N);
    title(['Q(x) eigenfunction series expansion with ', num2str(N), ' Terms']);
end
