function problem_2(terms_to_use)
    w = @(x) x.^3;
    figure; % Open a new figure window
    for index = 1:length(terms_to_use)
        lam = getlam(terms_to_use(index));
        Nx = 4*length(lam); if Nx<100, Nx=100; end
        dx = 4/Nx; xx = 0:dx:4;
        gg = QQ(xx); % Updated to use the new function QQ
        frep = zeros(size(xx));
        yn = zeros(size(xx));

        for ll = lam
            yn(1) = ll/2;
            yn(2:end) = besselj(1, ll*xx(2:end))./xx(2:end);
            NN = simp_int(yn.^2.*w(xx),dx);
            fhat = simp_int(yn.*gg.*w(xx),dx);
            frep = frep + fhat*yn/NN;
        end

        subplot(2, 2, index); % Create a subplot for each term
        plot(xx, frep, ':', xx, gg, '-k', 'linewidth', 2);
        title(sprintf('Using %d terms of eigens', terms_to_use(index)));
    end
end

function foo = simp_int(a, dx)
    foo = dx/3*(a(1)+a(end)+4*sum(a(2:2:end-1))+2*sum(a(3:2:end-2)));
end

function lams = getlam(N)
    lams = zeros(1, N); iroot = 0;
    ff = @(l) 0.5*(besselj(0, 4*l) - besselj(2, 4*l));
    x = 0; dx = 1/10;
    while iroot < N
        while ff(x) * ff(x+dx) > 0, x = x + dx; end
        iroot = iroot + 1;
        lams(iroot) = fzero(ff, x);
        x = x + 2*dx;
    end
end

function foo = QQ(x)
    % Nested function for QQ
    foo = zeros(size(x));
    foo(x >= 1 & x < 2) = abs(sin(6 * (x(x >= 1 & x < 2) - 1)));
    foo(x >= 2 & x < 3) = sqrt(x(x >= 2 & x < 3) - 2) - sin(6);
    foo(x >= 3) = 0;
end

%problem_2([5, 10, 20, 40]);
