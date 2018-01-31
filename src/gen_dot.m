%
% generate dot representation of undirected or directed graph
%
function gen_dot(A,F,Dir)

scale = 0.1;
[N,n] = size(A);

FN=fopen(F,'w');

if (Dir > 0),
	Op = '->';
	fprintf(FN,'digraph gr {\n');
	else
	Op = '--';
	fprintf(FN,'graph gr {\n');
end;

fprintf(FN,'node [shape=point, label="\\N"];\n');
%fprintf(FN,'graph [size="1,2", bb="0,0,30,100"];\n');
%fprintf(FN,'graph [size="1,2"];\n');
%fprintf(FN,'graph [size="10,20"];\n');
%fprintf(FN,'graph [size="10,20, weight=0.1, ranksep=5"];\n');
%fprintf(FN,'graph [size="10,20", weight=0.1, ranksep=5];\n');
%fprintf(FN,'graph [size="10,20", ranksep=5];\n');
%fprintf(FN,'graph [size="3,5", ranksep=5];\n');
%fprintf(FN,'size="3,5";\n');
%fprintf(FN,'ranksep=5;\n');
%fprintf(FN,'width=50;\n');
for i=1:N,
	fprintf(FN,'n_%d;\n', i);
%	fprintf(FN,'n_%d [width=%f];\n', i,scale*sum(A(i,:)));
end;
for i=1:N,
	if (Dir > 0),
		UL=N;
	else
		UL=i-1;
	end;
	for j=1:UL,
		if (A(i,j) > 0),
%			fprintf(FN,'n_%d %s n_%d;\n',i,Op,j);
			fprintf(FN,'n_%d %s n_%d[arrowsize=0.3];\n',i,Op,j);
		end;
	end;
end;
fprintf(FN,'}\n');
fclose(FN);

