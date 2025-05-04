function [U,UXT] = presentation_1

w = 10;
N = 40;
t = 2;

%Divide the interval into grid points
h = 1/N;
Nx = linspace(0,1,N+1);

%Peclet and mesh peclet number
MeshPe = (abs(w)*h)/2;

%Prep RHS vector f and solution vector u with boundary conditions
f = zeros(N-1,1);
u = zeros(N+1,1);

%Terms of the matrix given by centred approximations
a = 1/h^2 + w/(2*h);
b = 2/(h^2);
c = 1/(h^2) - w/(2*h);

%Form a tridiagonal matrix
A = spdiags([-a b -c],-1:1,N-1,N-1);

%Create time step and initial condition
tao = 0.2;
dt = h/t;
u_n = ones(N-1,1);

%Video maker
% v = VideoWriter("time_problems.avi");
% open(v)

%Colour gradient
% cmap = flip(copper(61),1);
% set(gca(),'ColorOrder',cmap)
% hold on


U = [];
UXT = [];

for n=0:1:(tao/dt)

    hold off
    
    %Plot approximate solution
    u(2:N) = u_n;
    plot(Nx,u,"*--",'MarkerSize',8, "Color", "r")
    U = [U ; transpose(u)];
    
    %Plot exact solution
    hold on
    uxt = exact_problem_II_static(w,N,t,n);
    plot(Nx, uxt, "Color", 'b', LineWidth=1)
    UXT = [UXT ; uxt];

    %Edit the graph
    axis([0 1 -1 1])
    grid on
    set(gca,'FontSize',18)
    title('Crank-Nicolson', ...
    ['w = ',num2str(w) ', h = ',num2str(h) ', t = ',num2str(dt*n)], 'FontSize',16);

    %Add frame to video
    % frame = getframe(gcf);
    % writeVideo(v,frame)

    pause(1)

    % if (0 <= n) && (n <= 1)
    %     %Update time step
    %     u_n1 = ((dt)^-1.*eye(N-1)+A)\(f+(dt)^-1.*u_n); %Implicit Euler Scheme
    %     u_n = u_n1;
    % else
    %     %Update time step
    %     DT = (dt)^-1.*eye(N-1);
    %     u_n1 = (DT+0.5.*A)\(f-(0.5.*A-DT)*u_n); %Crank-Nicolson scheme
    %     u_n = u_n1;
    % end

    DT = (dt)^-1.*eye(N-1);
    u_n1 = (DT+0.5.*A)\(f-(0.5.*A-DT)*u_n); %Crank-Nicolson scheme
    u_n = u_n1;

end


end

