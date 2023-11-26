function y=f_ex1_check()
% ======= bessel trick solution ============
xx = linspace(1,5,201);
u1 = @(x) besseli(0,2*sqrt(x)); % potential error
u2 = @(x) besselk(0,2*sqrt(x)); % potential error
CC = [u1(1) u2(1); u1(5) u2(5)]\[-10;-18];
besseltrick=@(x) CC(1)*u1(x) + CC(2)*u2(x);

plot(xx, besseltrick(xx), 'DisplayName', 'Bessel Trick'); hold on
% ======= series solution ============
xx = linspace(1,5,51);
u1 = @(x) f_ex1(x,1);
u2 = @(x) f_ex1(x,2);
CC = [u1(1) u2(1); u1(5) u2(5)]\[-10;-18];
seriessol=@(x) CC(1)*u1(x) + CC(2)*u2(x);
plot(xx, seriessol(xx), '+', 'DisplayName', 'Series Solutuion'); hold off
%--------------------------------------------------------------------------

function y=f_ex1(x,which)
% ======= (1) ==================
%ff =@(k) -1/((k-2)*(k-1)+(2*(k-1))); % most likely potential error
ff =@(k) -1/(k*(k+1))
xmax = max(abs(x));
NN=2;
while abs(xmax*ff(NN))>1e-3, NN=NN+1; end
y1 = ones(size(x));
for kk=NN:-1:1
    y1 = 1 + ff(kk)*x.*y1;
end
if which==1, y=y1; return; end
% === second solution ==========
c0 = 1; cc=zeros(1,NN); cc(1) = ff(1)*c0';
d0 = 1; dd=zeros(1,NN); dd(1) = (d0 - 2*cc(1)*1)/(1)^2; % potential error
for kk=2:NN
    cc(kk) = ff(kk)*cc(kk-1);
    dd(kk) = (-dd(kk-1) - 2*cc(kk)*kk) / (kk^2+kk); % potential error
end
dd = [dd(end:-1:1) d0];
y = log(x).*y1 + polyval(dd,x);

% Add title and axis labels
title('Comparison of Bessel Trick and Series Solutions');
xlabel('x');
ylabel('y');
legend();