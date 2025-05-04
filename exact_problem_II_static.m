function uxt = exact_problem_II_static(w,N,timestep,i)

h = 1/N;
dt = h/timestep;

x=linspace(0,1,N+1);

t = dt*i;    

Vsum = 0;

if t == 0
    uxt = ones(1,N+1);
else
    for n = 1:30000
        Bn = ((8*exp(-0.5*w)*pi*n)*(exp(0.5*w)-(-1)^n))/(w^2+(2*pi*n)^2)*sin(n*pi*x)*exp(-(n*pi)^2*t);
        Vsum = Vsum + Bn;
    end

    vp = exp(0.5*w*(x-0.5*w*t)).*Vsum;

    uxt = vp;
end
 
end



