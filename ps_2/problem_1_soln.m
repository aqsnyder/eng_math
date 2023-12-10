%===================================
% List of Nterms to run
NtermsList = [5, 20, 40];

% Execute problem_1 for each Nterms
for N = NtermsList
    figure;
    problem_1(N);
    title(['Q(x) eigenfunction series expansion with ', num2str(N), ' Terms'])
end