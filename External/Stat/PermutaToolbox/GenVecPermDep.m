function [vp1,vp2] = GenVecPermDep(N)

perm = randperm(2*N);
vp1 = perm(1:N);

[temp,ind] = setdiff(perm,vp1);
%vp2 = perm(sort(ind));
for i = 1:N,
    if vp1(i) <= N 
        vp2(i)=vp1(i)+N;
    else
        vp2(i)=vp1(i)-N;
    end;
end;
vp1=vp1';
vp2=vp2';

