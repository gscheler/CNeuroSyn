function [m_isi, s_isi] = calc_isi_off(sim,off,do_plot)

for i=1:sim.N_nn,
	is=find(sim.instrument.allvm(1,i,off:end)>sim.activity_thr);
	iisi = [is;sim.T_upd] - [0;is];
	if (do_plot),
		figure;
		hist(iisi,10:10:60);
		end;
	xx = find(iisi < sim.isi_max & iisi >sim.isi_min);
	if (isempty(xx)),
		m_isi(i) = sim.isi_max;
		s_isi(i) = sim.isi_max;
	else
		m_isi(i) = mean(iisi(xx));
		s_isi(i) = std(iisi(xx));
		end;
	end;
