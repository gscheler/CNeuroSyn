
1) NO direction

load the LN and R_hist

AA_R = A_R|A_R';

gen_dot(AA_LN(40:89,40:89),'gr_LN_hist.dot',0);
gen_dot(AA_R(131:180,131:180),'gr_R_hist.dot',0);`

with goal: roughly same # of edges in subgraph

twopi -Tpng gr_LN_hist.dot > gr_LG1_50.png
twopi -Tps gr_LN_hist.dot > gr_LG1_50.ps
twopi -Tps gr_R_hist.dot > gr_RG_50.ps
twopi -Tpng gr_R_hist.dot > gr_RG_50.png

directed graph:
2) with direction:

A_LN=A_LN & ~eye(210);
gen_dot(A_LN(40:89,40:89),'gr_LN_hist2.dot',1);
gen_dot(A_R(131:180,131:180),'gr_R_hist2.dot',1);
imagesc(A_R-A_R')
AA=A_R(131:180,131:180);
imagesc(AA-2*AA')
AA(35,1)
A_LN=A_LN & ~eye(210);
gen_dot(A_LN(40:89,40:89),'gr_LN_hist2.dot',1);
gen_dot(A_R(131:180,131:180),'gr_R_hist2.dot',1);
imagesc(A_R-A_R')
AA=A_R(131:180,131:180);
imagesc(AA-2*AA')
AA(35,1)
A(35,1)=0;
imagesc(AA-2*AA')
AA(32,28)=0;
imagesc(AA-2*AA')
gen_dot(AA,'gr_R_hist2.dot',1);

% alyout with twopi


