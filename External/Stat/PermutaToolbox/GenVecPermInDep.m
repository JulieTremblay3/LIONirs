function [vp1,vp2] = GenVecPermInDep(n1,n2)

perm = randperm(n1+n2);
vp1 = perm(1:n1);

[temp,ind] = setdiff(perm,vp1);
vp2 = perm(sort(ind));

vp1=vp1';
vp2=vp2';

