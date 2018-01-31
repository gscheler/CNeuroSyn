
function B = chop(A, mx)

od=sum(A);

keep=find(od <= mx);

B = A(keep,keep);


