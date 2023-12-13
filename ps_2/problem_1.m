function Qfunc(N)
%===================================
Q = @(x) (0.*(x>=0 & x<1) + ...
          abs(sin(6*(x - 1))).*(x>=1 & x<2) + ...
          (sqrt(x - 2) - sin(6)).*(x>=2 & x<3) + ...
          0.*(x>=3 & x<=4));
dx = 1/100;
x = 0:dx:4;
%===================================
lam = ((2*(1:N)-1)*pi)/8;
w = ones(size(x));
q = Q(x);

Qfunc = zeros(size(x));
for n=1:N
    yn = sin(lam(n)*x);
    Nn = simpson(yn.^2 .*w, dx);
    fn = simpson(q.*yn.*w, dx);
    Qfunc = Qfunc + fn*yn/Nn;
end

plot(x, q, '-k','linewidth',2); hold on;
plot(x, Qfunc, '--'); hold off;

function foo = simpson(Q,dx)
    foo = dx/3*(Q(1) + Q(end) + 4*sum(Q(2:2:end-1)) + 2*sum(Q(3:2:end-2)));
end

end


