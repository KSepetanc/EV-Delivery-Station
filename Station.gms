* Parameters loading file
$include input_data.gms

* Equation declaration macro (c stands for constraint)
$macro c(name) equation name;\
                name

* Variables declaration
free variables OF,q(t);
integer variables n_bat(t,n),n_aux(t,n,m),x_in(t,n),x_out(t,k,n,m);
nonnegative variables b_deg(t),y(t),z,omega(t);

* Parameter used to speed-up GAMS model creation.
parameter N_cond(t,k,n);
N_cond(t,k,m)$(N_out(t-ord(k),k,m)>0)=1;

* Macros for dollar conditions
$macro k_cond(n,m) ((ord(m)>=ord(n)-K_dis(n)) and (ord(m)<=ord(n)+K_ch(n)))
$macro x_in_sparse1 (ord(n)>=N_min)
$macro x_out_sparse1(n,m) (ord(m)<=ord(n) and N_out(t,k,m)>0)
$macro x_out_sparse2 (N_cond(t,k,m)>0)
$macro energy (ord(k)>card(t)-ord(t) and ord(m)>=ord(n) and N_out(t,k,n)>0)
* Macro for hours is defined in input_data.gms.

* Limit minimum SoE cluster for outgoing EVs.
x_out.fx(t,k,n,m)$(ord(m)<21 and N_out(t,k,n)>0 and ord(m)>=ord(n))=0;


* Model "S" - Station
c(S_profit)..                            OF =e= sum(t$hours,Lambda(t)*(q(t)+q(t+1)+q(t+2)+q(t+3)) ) + sum(t,b_deg(t))
                                              + z*Gamma + sum(t$hours,omega(t));

c(S_robust1(t))$hours..                  z + omega(t) =g= Delta_lambda(t)*y(t);

c(S_robust2_1(t))$hours..               -y(t) =l= q(t)+q(t+1)+q(t+2)+q(t+3);
c(S_robust2_2(t))$hours..                q(t)+q(t+1)+q(t+2)+q(t+3) =l= y(t);

* Constraints
c(S_clust1(t,n))..                       n_bat(t-1,n) + N_start(n)$(ord(t)=1) - sum((k,m)$x_out_sparse1(n,m),x_out(t,k,m,n)) + x_in(t,n)$x_in_sparse1 =e= sum(m$k_cond(n,m),n_aux(t,n,m));

c(S_clust2(t,n))..                       n_bat(t,n) =e= sum(m$k_cond(m,n),n_aux(t,m,n));

c(S_io1(t,k,n))$(N_out(t,k,n)>0)..       sum(m$(ord(m)>=ord(n)),x_out(t,k,n,m)) =e= N_out(t,k,n);

c(S_io2(t,n))$x_in_sparse1..             x_in(t,n) =e= sum((k,m)$x_out_sparse2,x_out(t-ord(k),k,m,m+(ord(n)-N_min)) );

c(S_end)..                               sum(n,n_bat(TEnd,n)*C_bat*(ord(n)-1)/(N_clus-1)) + sum((t,k,n,m)$energy,x_out(t,k,n,m)*C_bat*(ord(m)-ord(n)+N_min-1)/(N_clus-1))
                                                   =g= sum(n,N_start(n)*C_bat*(ord(n)-1)/(N_clus-1));

c(S_deg(t))..    b_deg(t) =e= sum(n,sum(m$(ord(m)<ord(n) and k_cond(n,m)),Delta*n_aux(t,n,m)*C_bat*(ord(n)-ord(m))/(N_clus-1)) );

c(S_q(t))..      q(t) =e= sum(n, sum(m$(ord(m)>ord(n) and k_cond(n,m)),n_aux(t,n,m)*C_bat/Eta(n,m) *(ord(m)-ord(n))/(N_clus-1))
                                -sum(m$(ord(m)<ord(n) and k_cond(n,m)),n_aux(t,n,m)*C_bat*Eta(n,m) *(ord(n)-ord(m))/(N_clus-1)) );

c(S_maxN(t))..   sum(n, sum(m$(ord(m)>ord(n) and k_cond(n,m)),n_aux(t,n,m))
                                +sum(m$(ord(m)<ord(n) and k_cond(n,m)),n_aux(t,n,m)) ) =l= N_ch;


model Station /All/;

option optcr=0.001;
option threads=4;
option mip = gurobi;

*Set variable bounds
parameter n_start_tot;
n_start_tot=sum(n,n_start(n));
n_bat.up(t,n)=n_start_tot;
n_aux.up(t,n,m)$k_cond(n,m)=n_start_tot;
x_in.up(t,n)=n_start_tot;
x_out.up(t,k,n,m)$x_out_sparse1(m,n)=n_start_tot;

* Use solver's barrier algorithm for more difficult models, e.g. a week or a month simulated.
file opt1 gurobi option file /gurobi.opt/;
put opt1;
put ''/;
*put 'method 2'/;
putclose;

Station.OptFile = 1;

solve Station using mip minimizing OF;

