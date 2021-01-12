* Parameters loading file
$include input_data.gms

* Equation declaration macro (c stands for constraint)
$macro c(name) equation name;\
                name

parameter qdam(t),qmpc(t),x_in_par(tl,n),qdelta(t),end_energy,x_in_par2(tl,n),N_start2(n),bdeg_mpc(t);
scalar time_now;
execute_load 'day_ahead_decisions.gdx' qdam=q.l;
x_in_par(tl,n)=0;
qmpc(t)=0;
bdeg_mpc(t)=0;
N_start2(n)=N_start(n);


* Variables declaration
free variables OF,q(t);
integer variables n_bat(t,n),n_aux(t,n,m),x_in(t,n),x_out(tl,k,n,m),b_deg(t);


* Parameter used to speed-up GAMS model creation.
parameter N_cond(tl,k,n);
N_cond(tl,k,m)$(N_out(tl-ord(k),k,m)>0)=1;

* Macros for dollar conditions
$macro k_cond(n,m) ((ord(m)>=ord(n)-K_dis(n)) and (ord(m)<=ord(n)+K_ch(n)))
$macro x_in_sparse1 (ord(n)>=N_min)
$macro x_out_sparse1(n,m) (ord(m)<=ord(n) and N_out(t,k,m)>0)
$macro x_out_sparse2 (N_cond(t,k,m)>0)
$macro energy (ord(k)>card(t)-ord(t) and ord(m)>=ord(n) and N_out(t,k,n)>0)
* Macro for hours is defined in input_data.gms.


* Limit minimum SoE cluster for outgoing EVs.
*x_out.fx(t,k,n,m)$(ord(m)<21 and N_out(t,k,n)>0 and ord(m)>=ord(n))=0;


* Model "S" - Station
c(S_max_energy)..                        OF =e=  sum(t$(ord(t)<time_now and hours(0)),IM_Lambda(t-3)*(qmpc(t-3)+qmpc(t-2)+qmpc(t-1)+qmpc(t)))
                                                +sum(t$(ord(t)=time_now and hours(1)),IM_Lambda(t)*(q(t)+q(t+1)+q(t+2)+q(t+3)))
                                                +sum(t$(ord(t)=time_now and hours(2)),IM_Lambda(t-1)*(qmpc(t-1)+q(t)+q(t+1)+q(t+2)))
                                                +sum(t$(ord(t)=time_now and hours(3)),IM_Lambda(t-2)*(qmpc(t-2)+qmpc(t-1)+q(t)+q(t+1)))
                                                +sum(t$(ord(t)=time_now and hours(0)),IM_Lambda(t-3)*(qmpc(t-3)+qmpc(t-2)+qmpc(t-1)+q(t)))
                                                +sum(t$(ord(t)>time_now and hours(1)),IM_Lambda(t)*(q(t)+q(t+1)+q(t+2)+q(t+3)))
                                                +sum(t$(ord(t)>=time_now),b_deg(t)) +sum(t$(ord(t)<time_now),bdeg_mpc(t)) ;


c(S_deg(t))$(ord(t)>=time_now)..         b_deg(t) =e= sum(n,sum(m$(ord(m)<ord(n) and k_cond(n,m)),Delta*n_aux(t,n,m)*C_bat*(ord(n)-ord(m))/(N_clus-1)) );

c(S_clust1(t,n))$(ord(t)>=time_now)..    n_bat(t-1,n)$(ord(t)>time_now) + N_start2(n)$(ord(t)=time_now) - sum((k,m)$x_out_sparse1(n,m),x_out(t,k,m,n)) + x_in(t,n)$x_in_sparse1 + x_in_par(t,n)=e= sum(m$k_cond(n,m),n_aux(t,n,m));

c(S_clust2(t,n))$(ord(t)>=time_now)..    n_bat(t,n) =e= sum(m$k_cond(m,n),n_aux(t,m,n));

c(S_io1(t,k,n))$(N_out(t,k,n)>0 and ord(t)>=time_now)..  sum(m$(ord(m)>=ord(n)),x_out(t,k,n,m)) =e= N_out(t,k,n);

c(S_io2(t,n))$(x_in_sparse1 and ord(t)>=time_now)..      x_in(t,n) =e= sum((k,m)$(x_out_sparse2 and (ord(t)-ord(k)>=time_now)),x_out(t-ord(k),k,m,m+(ord(n)-N_min)) );


c(S_q(t))$(ord(t)>=time_now)..      q(t)+qdam(t) =e= sum(n, sum(m$(ord(m)>ord(n) and k_cond(n,m)),n_aux(t,n,m)*C_bat/Eta(n,m) *(ord(m)-ord(n))/(N_clus-1))
                                   -sum(m$(ord(m)<ord(n) and k_cond(n,m)),n_aux(t,n,m)*C_bat*Eta(n,m) *(ord(n)-ord(m))/(N_clus-1)) );

c(S_end)..                          sum(n,n_bat(TEnd,n)*C_bat*(ord(n)-1)/(N_clus-1)) + sum((t,k,n,m)$(energy and ord(t)>=time_now),x_out(t,k,n,m)*C_bat*(ord(m)-ord(n)+N_min-1)/(N_clus-1))
                                  + sum((tl,n)$(ord(tl)>96),x_in_par(tl,n)*C_bat*(ord(n)-1)/(N_clus-1)) =g= sum(n,N_start(n)*C_bat*(ord(n)-1)/(N_clus-1));

*c(S_maxN(t))$(ord(t)>=time_now)..   sum(n, sum(m$(ord(m)>ord(n) and k_cond(n,m)),n_aux(t,n,m))
*                                   +sum(m$(ord(m)<ord(n) and k_cond(n,m)),n_aux(t,n,m)) ) =l= N_ch;


model MPC /All/;

option intVarUp=0;
option optcr=0.0;
option threads=4;
option mip = gurobi;

MPC.OptFile = 1;


for (time_now=1 to 96,
solve MPC using mip minimizing OF;

N_start2(n)= sum(t$(ord(t)=time_now),n_bat.l(t,n));

qmpc(t)$(ord(t)=time_now)=q.l(t);
bdeg_mpc(t)$(ord(t)=time_now)=b_deg.l(t);

x_in_par(tl,n)$(ord(tl)>time_now) = x_in_par(tl,n) + sum((k,m)$((N_cond(tl,k,m)>0) and (ord(tl)-ord(k)=time_now)),x_out.l(tl-ord(k),k,m,m+(ord(n)-N_min)) );

*Activate to make EVs arrive earlier.
x_in_par(tl,n)$(ord(tl)=time_now+1)=x_in_par(tl,n)+x_in_par(tl+3,n);
x_in_par(tl,n)$(ord(tl)=time_now+4)=0;

*Activate to make EVs arrive with higher SoE.
*x_in_par2(tl,n)$(ord(tl)=time_now+1)=x_in_par(tl,n-1);
*x_in_par(tl,n)$(ord(tl)=time_now+1)=x_in_par2(tl,n);

);
