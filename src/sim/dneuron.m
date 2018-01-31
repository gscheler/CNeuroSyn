%
% IZH neuron - simple model with Euler integration
%
function [v,u,ind] = dneuron(v,u,I_syn,a,b,c,d)

dt = 1;

dv = 0.04*v.^2 + 5*v + 140 - u - I_syn;

du = a.*(b.*v - u);

v = v + dv*dt;
u = u + du*dt;

ind=find(v>30);
v(ind) = c(ind);
u(ind) = u(ind) +d(ind);


