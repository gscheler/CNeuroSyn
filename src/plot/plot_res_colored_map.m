%
% plot 2D density x sigma* --> s
% must have loaded
%	res for the multiple graphs
%
	% number of grid-points
N_lin=10;
N_lin=20;

x=res.sigma_star';
y=res.kden';
z=res.sync_G';

%x=[res_R.sigma_star res.sigma_star]';
%y=[res_R.kden res.kden]';
%z=[res_R.sync_G res.sync_G]';

%xlin=linspace(min(x),max(x),N_lin);
%ylin=linspace(0.01,0.1,N_lin);
xlin=linspace(1.1,3,N_lin);
ylin=linspace(0.01,0.1,N_lin);

[X,Y] = meshgrid(xlin,ylin);
f=scatteredInterpolant(x,y,z);

Z=f(X,Y);

figure
mesh(X,Y,Z);

Z(Z<0)=0;
figure;
colormap jet
pcolor(X,Y,Z);
hold on
shading interp

plot(x,y,'k.')

axis([1.1,3,0.01,0.1])

xlabel('\sigma^*')
ylabel('density')

set(gca,'fontsize',12);
h=colorbar;
set(h,'XTick',[0:0.1:0.5]);
