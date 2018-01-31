% Plot raster plots for set of neurons
% 
% plot_syn_response.m
%	sim: = sim structure 
%	nn_inputs = NN inputs
%	LW = line width for plotting
%
function [isi_grp_mean, isi_grp_std] = plot_nfig11_2018(FN, sim, nn_inputs, LW, nn_pars, off)

[m_isi, s_isi] =calc_isi_plot(sim, off, 0);

N =sim.N_nn;

figure;

%--------------------------------------------
firact = zeros(1,N);

m=colormap('jet');
%m=colormap('cool');
m=zeros(64,3);
c=1;
%--------------------------------------------
ti=1:sim.T_upd-off+1;
RP=randperm(210);
hold off;
for ii=1:sim.N_nn,
	i = RP(ii);
	col = 1+4*floor(i/30);
        sp=find(sim.instrument.allvm(1,i,off:end) > sim.activity_thr);
	fireact(i)=length(sp);
        if (length(sp) > 0),
                spp=zeros(1,length(sp))+ii;
%                plot(ti(sp),spp,'k.','MarkerSize',8);
                plot(ti(sp)-1000,spp,'.','MarkerSize',4,'color',m(col,:));
		hold('on');
        end;
end;

%set(gca,'visible','off');


%plot([1800,1900],[-10,-10],'k','linewidth',5);
%text(1900,-10,'100ms');
%axis(gca,[999,2000,-10,210]);

set(gca,'Fontsize',12);
xlabel('ms');
ylabel('Neuron Number');

axis(gca,[999-1000,2000-1000,1,210]);




return;




s=1;
dl=30;
for i=1:7,
%	sigma_2 = std(m_isi(s+(i-1)*dl:s+i*dl-1));
%	mmm = mean(m_isi(s+(i-1)*dl:s+i*dl-1));

	msel = m_isi(s+(i-1)*dl:s+i*dl-1);
	ms = msel(find(msel < 500));
	mmm = mean(ms);
	sigma_2 = std(ms);

if(1==0),
	text_isi = sprintf('I=%.1f[%.1f]',mmm, sigma_2);
	text(-150,s+dl*(i-1)+10,text_isi, 'FontSize',[12]);
end;

	isi_grp_mean(i) = mean(m_isi(s+(i-1)*dl:s+i*dl-1));
	isi_grp_std(i) = sigma_2;

if (1==0),
	n_fire = sprintf('%3d',...
	   mean(neuron_activity(sim.instrument.allvm(1,s+(i-1)*dl+5,:),off,sim)));
else
	na = fireact(s+(i-1)*dl:s+i*dl-1);
	n_fire = sprintf('%.1f [%.1f]', mean(na), std(na));
end;
if (1==0),
	text(sim.T_upd-off,s+dl*(i-1)+10,n_fire, 'FontSize',[12]);
end;
	end;

hold off;
axis([1,sim.T_upd-off,1,sim.N_nn+1]);
%set(gca,'Visible','off');

%subplot(3,1,3);
%sel=3;
%        sp=find(sim.instrument.allvm(1,sel,off:end)> sim.activity_thr);
%        if (length(sp) > 0),
%                spp=zeros(1,length(sp))+3;
%                plot(ti(sp),spp,'r.','MarkerSize',8);
%        end;

%--------------------------------------------

%fn_eps =sprintf('%s.eps', FN);
%print('-deps', fn_eps);
%fn_jpg =sprintf('%s.jpg', FN);
%print('-djpeg', fn_jpg);
