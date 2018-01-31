%
% calculate neuron-neuron spike correlation
%
function [m_sc, w_act] = calc_spikecorr2(sim, off, c_max,W)

m_sc = zeros(sim.N_nn,sim.N_nn);
s_sc = zeros(sim.N_nn,sim.N_nn);
w_act = zeros(sim.N_nn, sim.T_upd);
act = zeros(1, sim.N_nn);

%----------------------------------------------------------

for nn=1:sim.N_nn,
  for i=off:W:sim.T_upd-W,
	a = find(sim.instrument.allvm(1,nn,i:i+W-1)>sim.activity_thr);
	if (a>0),
		w_act(nn,i) = 1;
		act(nn) = act(nn)+1;
	end;
  end;
end;
	
for nn=1:sim.N_nn,
  for i=find(w_act(nn,:)>0),
	lm=max(1,i-c_max);
	rm=min(sim.T_upd-W,i+c_max);
	for j=nn+1:sim.N_nn,
		a=find(w_act(j,lm:rm)>0);
		if (a>0),
			m_sc(nn,j) = m_sc(nn,j) + 1;
			m_sc(j,nn) = m_sc(j,nn) + 1;
		end;
	end;
  end;
end;

