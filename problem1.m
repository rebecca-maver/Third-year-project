function u = problem1(w,N)

%Divide the interval into grid points
h = 1/N;
Nx = linspace(0,1,N+1);

%Peclet and mesh peclet number
MeshPe = (abs(w)*h)/2;

%Prep a matrix A, RHS vector f and boundary conditions
f = h^2.*ones(N-1,1);
u = zeros(N+1,1);

%Terms of the matrix given by centred approximations
a = 1 + (w*h)/2;
b = 2;
c = 1 - (w*h)/2;

%Form the tridiagonal matrix
% tic
% A = zeros(N-1);
% A(1,1) = b;
% for i=2:N-1
%     A(i,i-1) = -a;
%     A(i,i) = b;
%     A(i-1,i) = -c;
% end
% toc

tic
A = spdiags([-a b -c],-1:1,N-1,N-1);
toc

%Solve for a vector of approximates u
u(2:N) = A\f;

%Plot for x in [0, 1]
ux = @(x) (1/w)*(x - (1-exp(w*x))/(1-exp(w)));
fplot(ux, [0 1], "k", 'LineWidth',1)
hold on
plot(Nx,u,"*:",'MarkerSize',10,'LineWidth',2)
grid on
set(gca,'FontSize',18)
title(['Mesh Peclet number : ' num2str(MeshPe)], ...
    ['w = ',num2str(w) ', h = ',num2str(h)], 'FontSize',22);


end

