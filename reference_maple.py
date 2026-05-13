#http://www.mathematik.uni-kl.de/~numberfieldtables/
read "classnmb":
with(numtheory):
read "SF2.4v.txt":
with(SF):


words(0):kernelopts(printbytes=false):


#This function computes the Hurwitz-Kronecker numbers (using the definition that takes units of the orders into account).
#Note that this function gives an error if dscr is negative and 0 or 3 mod 4.
HKclass:=proc(dscr) local fac_dscr,cnd,fund_dscr,i,j,k,res,div_cnd,clnmb,fac_subcnd:
fac_dscr:=ifactors(dscr)[2];
#
cnd:=1: fund_dscr:=-1:
for i from 1 to nops(fac_dscr) do
if fac_dscr[i,2] mod 2=0 then
cnd:=cnd*fac_dscr[i,1]^(fac_dscr[i,2]/2):
else
cnd:=cnd*fac_dscr[i,1]^((fac_dscr[i,2]-1)/2):
fund_dscr:=fund_dscr*fac_dscr[i,1]:
fi:
od:
#
if fund_dscr mod 4=2 or fund_dscr mod 4=3 then
fund_dscr:=4*fund_dscr: cnd:=cnd/2:
fi:
#
res:=0:
div_cnd:=[op(divisors(cnd))];
for i from 1 to nops(div_cnd) do
clnmb:=LstClnmb[fund_dscr]*div_cnd[i]; fac_subcnd:=ifactors(div_cnd[i])[2]; k:=1:
#
if nops(fac_subcnd)>0 then
#
if fac_subcnd[1,1]=2 then k:=2:
if fund_dscr mod 8 =1 then
clnmb:=clnmb*(1-1/2):
else if fund_dscr mod 8=5 then
clnmb:=clnmb*(1+1/2):
fi:fi:fi:
#
for j from k to nops(fac_subcnd) do
clnmb:=clnmb*(1-legendre(fund_dscr,fac_subcnd[j,1])/fac_subcnd[j,1]):
od:#special treatment dscr=-1, dscr=-3
fi:
res:=res+clnmb: #print(div_cnd[i],clnmb):
od:
#
if fund_dscr=-3 then res:=res/6; else
if fund_dscr=-4 then res:=res/4; else
res:=res/2;
fi:fi:
res:
end:

#Small note:
#Consider the case the curves with 4 automorphisms and with r odd.  
#If p=3 mod 4, all 4 have beta=0.
#If p=1 mod 4, then take a,b such that (a+b*I)*(a-b*I)=p, and put pi=a+b*I, then the four curves have beta=pi+bar(pi), -pi-bar(pi), I*pi-I*bar(pi), -I*pi+I*bar(pi). 
#Reference page 307 of Ireland and Rosen.
#Consider the case the curves with 6 automorphisms and with r odd.  And p not 2. 
#If p=2 mod 3, all six have beta=0.
#If p=1 mod 3, then take a,b such that (a+b*w)*(a+b*w^2)=p, and put pi=a+b*w, then the six curves have beta=pi+bar(pi), -(pi+bar(pi)), w*pi+w^2*bar(pi), -(w*pi+w^2*bar(pi)), w^2*pi+w*bar(pi), -(w^2*pi+w*bar(pi)). 
#Reference, page2 304-307 of Ireland and Rosen.

print("count_A1q(p,r)");
lprint("This function counts, for each k-isogeny class, the number of k-isomorphism classes of elliptic curves defined over k, weighted by the reciprocal of their number of k-automorphisms. It stores the result in the variable cat(A1q,p^r), where p^r is the number of elements of k.");
#
count_A1q:=proc(p,r) local q,Res,amax,tst,a,dscr,res:
q:=p^r:
if q<10^7/4 then
#
Res:=table():
amax:=floor(2*sqrt(q));
tst:=0:
#
for a from 0 to amax do
#
if not(a mod p=0) then
dscr:=a^2-4*q: res:=HKclass(dscr): print([a],res,case1):
else
if (a=0 and r mod 2=1) then
dscr:=a^2-4*p: res:=HKclass(dscr): print([a],res,case2):
else
if (a^2=2*q and p=2 and r mod 2=1) then
res:=1/4: print([a],res,case3):
else
if (a^2=3*q and p=3 and r mod 2=1) then
res:=1/6: print([a],res,case4):
else
if (a^2=4*q and r mod 2=0 and p=2) then
res:=1/24: print([a],res,case5):
else
if (a^2=4*q and r mod 2=0 and p=3) then
res:=1/12: print([a],res,case6):
else
if (a^2=4*q and r mod 2=0) then
res:=1/24*(p+6-4*legendre(-3,p)-3*legendre(-4,p)):
if not(p mod 3=1) then res:=res-1/2+1/6: fi: #E with j=0 appears here (check).
if not(p mod 4=1) then res:=res-1/2+1/4: fi: #E with j=0 appears here (check).
print([a],res,case7):
else
if (a^2=q and r mod 2=0) then
res:=(1-legendre(-3,p))/6: print([a],res,case8):
else
if (a=0 and r mod 2=0) then
res:=(1-legendre(-4,p))/4: print([a],res,case9):
else
res:=0: print([a],0,emptycase):
fi:fi:fi:fi:fi:fi:fi:fi:fi:
#
Res[[a]]:=res:
Res[[-a]]:=Res[[a]]:
#
if a=0 then tst:=tst+Res[[a]]: else tst:=tst+2*Res[[a]]: fi:
##if not (Res[[a]]-eval(cat(A1q,q))[[a]])=0 then print("mistake"): fi:
od:
if not(tst=q) then print("Mistake!"): fi:
cat(A1q,q):=op(Res):
else print("q too large"): fi:
end:

#This function applies count_A1q to all k whose number of elements lies between "strt" and "fin".
count_A1q_all:=proc(strt,fin) local n,nn,p,r:
for n from strt to fin do
nn:=ifactors(n)[2]:
if nops(nn)=1 then
p:=nn[1,1]: r:=nn[1,2]: print("q"=p^r):
count_A1q(p,r):
fi:od:
end:

#This function applies count_A1q to all k whose number of elements lies in Q.
count_A1q_Q:=proc(Q) local n,nn,p,r:
for n from 1 to nops(Q) do
nn:=ifactors(Q[n])[2]:
if nops(nn)=1 then
p:=nn[1,1]: r:=nn[1,2]: print("q"=p^r):
cat(A1q,Q[n]):=count_A1q(p,r):
fi:od:
end:

print(); print("sk(k)");
print("This function computes the dimension of the space S_k of elliptic modular forms of weight k.");
sk:=proc(k)
local s: s:=0:
if k<0 then s:=0 elif irem(k,2)=1 then s:=0 elif k=2 then s:=-1
elif irem(k,12)=2 then
s:=floor(k/12)-1 else s:=floor(k/12) fi:
s:
end:

print();print("count_Sk(p,r,kk)");
lprint("This function computes the trace of Frobenius F_{p^r} on the Galois representation corresponding to the space S_k of elliptic modular forms of weight kk.");
count_Sk:=proc(p,r,kk) local pol,k,q,Res,amax,tst,a,dscr,res:
q:=p^r: k:=kk-2: pol:=foldl(`+`,0,seq(e1^(k-2*i)*(-1)^i*q^i*binomial(k-i,i),i=0..k/2));
if q<10^7/4 then
if kk mod 2=1 then res:=0: else
#
Res:=0:
amax:=floor(2*sqrt(q));
tst:=0:
#
for a from 0 to amax do
#
if not(a mod p=0) then
dscr:=a^2-4*q: res:=HKclass(dscr)*2*subs(e1=a,pol): #print([a],res,case1):
else
if (a=0 and r mod 2=1) then
dscr:=a^2-4*p: res:=HKclass(dscr)*subs(e1=a,pol): #print([a],res,case2):
else
if (a^2=2*q and p=2 and r mod 2=1) then
res:=1/4*2*subs(e1=a,pol): #print([a],res,case3):
else
if (a^2=3*q and p=3 and r mod 2=1) then
res:=1/6*2*subs(e1=a,pol): #print([a],res,case4):
else
if (a^2=4*q and r mod 2=0 and p=2) then
res:=1/24*2*subs(e1=a,pol): #print([a],res,case5):
else
if (a^2=4*q and r mod 2=0 and p=3) then
res:=1/12*2*subs(e1=a,pol): #print([a],res,case6):
else
if (a^2=4*q and r mod 2=0) then
res:=1/24*(p+6-4*legendre(-3,p)-3*legendre(-4,p)):
if not(p mod 3=1) then res:=res-1/2+1/6: fi: #E with j=0 appears here (check).
if not(p mod 4=1) then res:=res-1/2+1/4: fi: #E with j=0 appears here (check).
res:=res*2*subs(e1=a,pol): #print([a],res,case7):
else
if (a^2=q and r mod 2=0) then
res:=(1-legendre(-3,p))/6*2*subs(e1=a,pol): #print([a],res,case8):
else
if (a=0 and r mod 2=0) then
res:=(1-legendre(-4,p))/4*subs(e1=a,pol): #print([a],res,case9):
else
res:=0: #print([a],0,emptycase):
fi:fi:fi:fi:fi:fi:fi:fi:fi:
#
Res:=Res+res:
od: fi: else print("q too large"): fi:
-(Res+1):
end:

print(); print("count_Sk_all(strt,fin,k)");
lprint("This function applies count_Sk to all finite fields whose number of elements lies between strt and fin and returns a list of the traces.");
count_Sk_all:=proc(strt,fin,k) local n,nn,p,r,Res,Q:
Res:=[]: Q:=[]:
for n from strt to fin do
nn:=ifactors(n)[2]:
if nops(nn)=1 then
p:=nn[1,1]: r:=nn[1,2]: #print("q"=p^r):
Q:=[op(Q),p^r]:
Res:=[op(Res),count_Sk(p,r,k)]:
fi:od:
print("Q",Q);
Res:
end:

print();print("charpol_Sk(p,kk)");
lprint("This function computes the characteristic polynomial of F_p acting on S[kk].");
charpol_Sk:=proc(pp,k) local dim,L,E,f,qq:
dim:=sk(k): qq:=pp^(k-1):
L:=[seq(count_Sk(pp,i,k),i=1..dim)]: #print(L):
E:=subs([seq(eval(cat(p,i))=eval(L[i]),i=1..nops(L))],[seq(top(cat(e,j)),j=1..nops(L))]); print(E);
#
f:=foldl(`+`,0,x^(2*dim),seq(x^(2*dim-i)*(-1)^i*E[i],i=1..dim),seq(x^i*qq^(dim-i)*(-1)^i*E[i],i=1..(dim-1)),qq^dim);
end:


print();print("count_ak(p,r,kk)");
lprint("This function computes the trace of Frobenius F_{p^r} on the Galois representation corresponding to the space S_k of elliptic modular forms of weight kk.");
count_ak:=proc(p,r,kk) local pol,k,q,Res,amax,tst,a,dscr,res:
q:=p^r: 
if q<10^7/4 then
if kk mod 2=1 then res:=0: else
#
Res:=0:
amax:=floor(2*sqrt(q));
tst:=0:
#
for a from 0 to amax do
#
if not(a mod p=0) then
dscr:=a^2-4*q: res:=HKclass(dscr)*2*a^kk: #print([a],res,case1):
else
if (a=0 and r mod 2=1) then
dscr:=a^2-4*p: res:=HKclass(dscr)*a^kk: #print([a],res,case2):
else
if (a^2=2*q and p=2 and r mod 2=1) then
res:=1/4*2*a^kk: #print([a],res,case3):
else
if (a^2=3*q and p=3 and r mod 2=1) then
res:=1/6*2*a^kk: #print([a],res,case4):
else
if (a^2=4*q and r mod 2=0 and p=2) then
res:=1/24*2*a^kk: #print([a],res,case5):
else
if (a^2=4*q and r mod 2=0 and p=3) then
res:=1/12*2*a^kk: #print([a],res,case6):
else
if (a^2=4*q and r mod 2=0) then
res:=1/24*(p+6-4*legendre(-3,p)-3*legendre(-4,p)):
if not(p mod 3=1) then res:=res-1/2+1/6: fi: #E with j=0 appears here (check).
if not(p mod 4=1) then res:=res-1/2+1/4: fi: #E with j=0 appears here (check).
res:=res*2*a^kk: #print([a],res,case7):
else
if (a^2=q and r mod 2=0) then
res:=(1-legendre(-3,p))/6*2*a^kk: #print([a],res,case8):
else
if (a=0 and r mod 2=0) then
res:=(1-legendre(-4,p))/4*a^kk: #print([a],res,case9):
else
res:=0: #print([a],0,emptycase):
fi:fi:fi:fi:fi:fi:fi:fi:fi:
#
Res:=Res+res:
od: fi: else print("q too large"): fi:
Res:
end:


###WORKING!!!
bincoeff_lucas:=proc(m,n,ell) local N,M,res:
if n>m then res:=0: else 
M:=convert(m,base,ell);
N:=convert(n,base,ell);
N:=[op(N),seq(0,i=1..nops(M)-nops(N))]:
res:=foldl(`*`,1,seq(binomial(M[i],N[i]) mod ell,i=1..nops(M))) mod ell:
fi:
res:
end:

#lprint("This function computes the trace modulo m of Frobenius F_{p^r} on the Galois representation corresponding to the space S_k of elliptic modular forms of weight kk.");
count_Sk_mod:=proc(p,r,kk,m) local Pol,TrE,k,q,Res,amax,tst,a,dscr,res:
q:=p^r: k:=kk-2:
if isprime(m) then Pol:=[seq((-1)^i*q&^i*bincoeff_lucas(k-i,i,m) mod m,i=0..k/2)];
else Pol:=[seq((-1)^i*q&^i*binomial(k-i,i) mod m,i=0..k/2)]; fi: #print(Pol);
#
if q<10^7/4 then
if kk mod 2=1 then res:=0: else
#
Res:=0:
amax:=floor(2*sqrt(q));
tst:=0:
#
for a from 0 to amax do
#
TrE:=foldl(`+`,0,seq(Pol[i+1]*a&^(k-2*i) mod m,i=0..k/2-1))+Pol[k/2+1] mod m:
#
if not(a mod p=0) then
dscr:=a^2-4*q: res:=modp(HKclass(dscr)*2*TrE,m): #print([a],res,case1):
else
if (a=0 and r mod 2=1) then
dscr:=a^2-4*p: res:=modp(HKclass(dscr)*TrE,m): #print([a],res,case2):
else
if (a^2=2*q and p=2 and r mod 2=1) then
res:=modp(1/4*2*TrE,m): #print([a],res,case3):
else
if (a^2=3*q and p=3 and r mod 2=1) then
res:=modp(1/6*2*TrE,m): #print([a],res,case4):
else
if (a^2=4*q and r mod 2=0 and p=2) then
res:=modp(1/24*2*TrE,m): #print([a],res,case5):
else
if (a^2=4*q and r mod 2=0 and p=3) then
res:=modp(1/12*2*TrE,m): #print([a],res,case6):
else
if (a^2=4*q and r mod 2=0) then
res:=modp(1/24*(p+6-4*legendre(-3,p)-3*legendre(-4,p)),m):
if not(p mod 3=1) then res:=modp(res-1/2+1/6,m): fi: #E with j=0 appears here (check).
if not(p mod 4=1) then res:=modp(res-1/2+1/4,m): fi: #E with j=0 appears here (check).
res:=modp(res*2*TrE,m): #print([a],res,case7):
else
if (a^2=q and r mod 2=0) then
res:=modp((1-legendre(-3,p))/6*2*TrE,m): #print([a],res,case8):
else
if (a=0 and r mod 2=0) then
res:=modp((1-legendre(-4,p))/4*TrE,m): #print([a],res,case9):
else
res:=0: #print([a],0,emptycase):
fi:fi:fi:fi:fi:fi:fi:fi:fi:
#
Res:=modp(Res+res,m):
od: fi: else print("q too large"): fi:
modp(-(Res+1),m):
end: