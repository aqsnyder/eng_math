function problem_2(Nterms)
    % Main function

    nint = 3 * max(Nterms);
    dx = 4 / nint; 
    xx = 0:dx:4;
    lam = get_roots(max(Nterms));

    for icase = 1:length(Nterms)
        NN = Nterms(icase);
        disp(NN);
        qq = zeros(size(xx));
        for nn = 1:NN
            YY = besselj(0, lam(nn) * xx);
            ww = xx;
            Nhat = simpson(ww .* YY .* YY, dx);
            Qhat = simpson(ww .* YY .* QQ(xx), dx);
            qq = qq + Qhat / Nhat * YY;
        end
        subplot(2, 2, icase);
        plot(xx, qq); 
        hold on;
        plot(xx, QQ(xx), 'k-', 'linewidth', 2); 
        hold off;
        title(sprintf('# terms = %d', NN));
    end

    function foo = QQ(x)
        % Nested function for QQ
        foo = zeros(size(x));
        foo(x >= 1 & x < 2) = abs(sin(6 * (x(x >= 1 & x < 2) - 1)));
        foo(x >= 2 & x < 3) = sqrt(x(x >= 2 & x < 3) - 2) - sin(6);
        foo(x >= 3) = 0;
    end

    function l = simpson(f, dx)
        % Nested function for Simpson's rule
        NN = length(f);
        l = dx / 3 * (f(1) + f(NN) + 4 * sum(f(2:2:NN - 1)) + 2 * sum(f(3:2:NN - 2)));
    end

    function lams = get_roots(NN)
        % Nested function for finding roots
        ff = @(x) besselj(0, 4 * x);
        dx = 1 / 10; x = 0;
        lams = zeros(1, NN);
        iroot = 0;
        while iroot < NN
            while ff(x) * ff(x + dx) > 0, x = x + dx; end
            iroot = iroot + 1;
            lams(iroot) = fzero(ff, x);
            x = lams(iroot) + dx;
        end
    end
end

%problem_2([5, 10, 20, 40])
