$eolCom //
set t 15 min time periods /t1*t96/;      //Time period can be increased up to t672, but increase the N_start so the model is feasable. (96)
set n SoE cluster /n1*n21/;
set w week price scenarios /w1*w52/;

alias(t,k);
alias(n,m);

singleton set TEnd(t);                   //Last time period set
TEnd(t)$(ord(t)=card(t))=t(t);

$macro hours (mod(ord(t),4)=1)

parameter N_start(n),N_out(t,k,n),Eta(n,m),Lambda_week(w,t),K_ch(n),K_dis(n),N_clus;
* Load data from GDX. GDX is created using xlsx_to_gdx.gms from excel.
$gdxin data.gdx
$load N_start,N_out,Eta,Lambda_week
$gdxin


* Maximum (dis)charging speed for each cluster.
K_ch(n)=sum(m$(ord(m)>ord(n) and Eta(n,m)),1);
K_dis(n)=sum(m$(ord(m)<ord(n) and Eta(n,m)),1);

N_clus=card(n);             //Number of clusters

parameter Lambda(t);        //Average electricity price, €/MWh
Lambda(t)$hours=sum(w,Lambda_week(w,t))/card(w);

parameter Gamma /8/;        //Uncertanty budget
parameter Delta_lambda(t);  //Price change interval, €/MWh
Delta_lambda(t)$hours=smax(w,abs(Lambda_week(w,t)-Lambda(t)));

parameter C_bat/0.044/;     //Battery capacity, MWh
parameter Delta/60/;        //Battery degradation cost, ï¿½/MWh
parameter N_min/5/;         //Vheicles minimum returning SoE cluster, 5th cluster -> 20% SoE
parameter N_ch/200/;         //Maximum number of vehicles that charge simultaneously.

* Shift N_out to match paper definition.
         //Offset time for handling and uncertinty
         //Change consumption to SoE requirements
         //Maximum data consumption cluster + N_min shouldn't be larger than N_clust+1
N_out(t,k,n)=N_out(t,k-3,n-(N_min-1));

* Scale the starting number of EVs.
n_start(n)=round(1*n_start(n));

