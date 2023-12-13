NtermsList = [5, 20, 40];

numSubplots = length(NtermsList);
figure;

for i = 1:numSubplots
    N = NtermsList(i);
    subplot(1, numSubplots, i);
    problem_1(N);
    title(['Q(x) eigenfunction series expansion with ', num2str(N), ' Terms']);
end
