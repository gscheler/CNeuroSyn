%
%  plot a Gaussian curves for the various groups of neurons
% given by s(),mi(), ma()
%
%
function plot_gaussian_ax(fn, m, s,mi,ma)


xmin=min(m)-max(s);
xmax=max(m)+max(s);
%x=-2:0.005:2;
x=xmin:(xmax-xmin)/500:xmax;

figure;
cmap=colormap('jet');
%cmap=colormap('winter');
%cmap=colormap('cool');
cmap=zeros(64,3);
%cmap(:,3) = linspace(1,0,64);
%c=10;
c=1;

for i=1:length(m),
	y=(1/sqrt(2*pi*s(i)))*exp(-((x-m(i)).^2)/s(i));
%	plot(x,y,'*');
	plot(x,y,'k','Linewidth',2,'color',cmap(c,:));
	hold on;
	fprintf('GROUP %d with mean: %f Hz\n',i,1000/m(i))
	c=c+4;
%	c=c+7;
	end;

%set(gca,'XTick',[50,100,143,200,250]);
%set(gca,'XTickLabel',{20,10,7,5,4});
%set(gca,'XTick',[50,100,200,250]);
%set(gca,'XTickLabel',{'20Hz/50ms','10Hz/100ms','5Hz/200ms','4Hz/250ms'});

xx=[50,100,150,200,250,300];
set(gca,'XTick',xx);
set(gca,'XTickLabel',round(1000./xx));

	%
	% neurons (30 Neurons per group)
	%
%set(gca,'YTick',[0.033,0.1,0.1667]);
%set(gca,'YTickLabel',{1,3,5});

ylabel('P','rotation',0);
%ylabel('number of neurons');
xlabel('Hz');
%xlabel('ms');
set(gca,'Fontsize',12);

axis([mi,ma,0,0.18]);

%fn_eps = sprintf('%s.eps', fn);
%print('-depsc',fn_eps);
%fn_png = sprintf('%s.png', fn);
%print('-dpng',fn_png);

