%
%	produce the degree histograms for the LG and RG graph
%
load LN_hist.mat
load RG_hist.mat


ID_LG = sum(A_LN');
ID_RG = sum(A_R');
OD_LG = sum(A_LN);
OD_RG = sum(A_R);

%ID_LG=ID_LG(ID_LG>0);
%OD_LG=OD_LG(OD_LG>0);
%ID_RG=ID_RG(ID_RG>0);
%OD_RG=OD_RG(OD_RG>0);

md=max([ID_LG,ID_RG,OD_LG,OD_RG]);


figure;

%[N,X]=hist([(OD_RG+ID_RG)',(OD_LG+ID_LG)'],1:md);
[N,X]=hist([(OD_RG)',(OD_LG+ID_LG-7)'],1:md);
h=bar(X,N,'barwidth',1,'edgecolor','none');
%hist([OD_RG',OD_LG'],30); %1:md);
%axis([0,60,0,60])



foo
[N,X]=hist(OD_RG',1:md);
h=bar(X,N,'barwidth',1,'edgecolor','none','facecolor','red','facealpha',0.5);
hold on
[N,X]=hist(OD_LG',1:md);
h=bar(X,N,'barwidth',1,'edgecolor','none','facecolor','blue','facealpha',0.5);

legend({'RG', 'LG1'})
xlabel('out-degree');
ylabel('number of connections');
set(gca,'fontsize',14');

%saveas(gcf,'hist-outdegree-RG-LG1','fig');
%saveas(gcf,'hist-outdegree-RG-LG1','png');

%axis([0,md,0,40])

figure;
%md=max([ID_LG,ID_RG]);
[N,X]=hist(ID_RG',1:md);
h=bar(X,N,'barwidth',1,'edgecolor','none','facecolor','red','facealpha',0.5);
hold on
[N,X]=hist(ID_LG',1:md);
h=bar(X,N,'barwidth',1,'edgecolor','none','facecolor','blue','facealpha',0.5);

legend({'RG', 'LG1'})
xlabel('in-degree');
ylabel('number of connections');
set(gca,'fontsize',14');
%axis([0,md,0,40])

%saveas(gcf,'hist-indegree-RG-LG1','png');
%saveas(gcf,'hist-indegree-RG-LG1','fig');

