* BIHAR PMP Agricultural Sector Risk Model  [B-PMP-ASRM] v01: May 2024
* Maxwell Mkondiwa
* Three similution
* (1)  Baseline model
* (2)  Calibration model
* (3)  simulation of different level of volatility performed by inflating the baseline V-COV matrix by an exogeous constant


* CARA, no RTS, farm specific alpha

SET n farms /1*70/;

SET j activities /Gram                   Gram
                  Jute                   Jute
                  Lentil                 Lentil
                  Maize                  Maize
                  Moong                  Moong
                  Potato                 Potato
                  Rapeseed_Mustard       Rapeseed_Mustard
                  Rice                   Rice
                  Sugarcane              Sugarcane
                  Tur                    Tur
                  Wheat                  Wheat /;

SET i resource /land/;

SET t years /t1*t12/;

Parameters  ha(n,j)    hectares
            prodt(n,j) production
            p(n,j)     output price
            ct(n,j)    costs per activity
            sp(t,j)    output price series;

$iftheni not exist "Bihar_PMP_Risk_Data.gdx"

$CALL GDXXRW.EXE Bihar_PMP_Risk_Data.xls par=ha rng=ha_avg!c1:n71 Rdim=1 Cdim=1

$GDXIN Bihar_PMP_Risk_Data.gdx
$LOAD ha
$GDXIN

Display ha;


$CALL GDXXRW.EXE Bihar_PMP_Risk_Data.xls par=prodt rng=prod_avg!c1:n71 Rdim=1 Cdim=1

$GDXIN Bihar_PMP_Risk_Data.gdx
$LOAD prodt
$GDXIN

Display prodt;


$CALL GDXXRW.EXE Bihar_PMP_Risk_Data.xls par=p rng=price!c1:n71 Rdim=1 Cdim=1

$GDXIN Bihar_PMP_Risk_Data.gdx
$LOAD p
$GDXIN

Display p;


$CALL GDXXRW.EXE Bihar_PMP_Risk_Data.xls par=ct rng=cost!c1:n71 Rdim=1 Cdim=1

$GDXIN Bihar_PMP_Risk_Data.gdx
$LOAD ct
$GDXIN

Display ct;


$CALL GDXXRW.EXE Bihar_PMP_Risk_Data.xls par=sp rng=rev_series!b1:m13 Rdim=1 Cdim=1

$GDXIN Bihar_PMP_Risk_Data.gdx
$LOAD sp
$GDXIN

Display sp;

    execute_unload "C:\Users\Maxwell\OneDrive\CIMMYT\Papers\PMPforEasternIndia\Data\PMP Risk Model\Bihar_PMP_Risk_Data.gdx" ha, prodt, p, ct, sp;
* execute_unload "D:\OneDrive\CIMMYT\Papers\PMPforEasternIndia\Data\GAMS_Model_Data\Bihar_PMP_Risk_Data.gdx" ha, prodt, p, ct, sp;

$else

    execute_load    "C:\Users\Maxwell\OneDrive\CIMMYT\Papers\PMPforEasternIndia\Data\PMP Risk Model\Bihar_PMP_Risk_Data.gdx" ha, prodt, p, ct, sp;
*execute_load    "D:\OneDrive\CIMMYT\Papers\PMPforEasternIndia\Data\GAMS_Model_Data\Bihar_PMP_Risk_Data.gdx" ha, prodt, p, ct, sp;
$endif

Parameter epsilon(n);
epsilon(n)=0.00000001;

Parameter c(n,j) cost per unit ;
c(n,j)$(prodt(n,j) gt 0)=(ct(n,j)*ha(n,j))/prodt(n,j);
display c;

Parameter r(n,j) yields;
r(n,j)$(ha(n,j) gt 0)=prodt(n,j)$(prodt(n,j) gt 0)/ha(n,j);
display r;

sp(t,j)=sp(t,j);
display sp;

*scaling
r(n,j)   = r(n,j);
p(n,j)   = p(n,j)*1000;
sp(t,j)  = sp(t,j)/1000;
c(n,j)   = c(n,j)/1000;

Parameter xbar(n,j) output level;
xbar(n,j)$(ha(n,j) gt 0)=ha(n,j)*r(n,j);
* Sugarcane is grown throughout the year
xbar(n,"Sugarcane")=xbar(n,"Sugarcane")/2;
display xbar;


Parameter msp(j);
msp(j)=sum(t, sp(t,j))/sum(t $ sp(t,j),1);
display msp;

ALIAS(j,k,jj,kk);

Parameter S(j,k) variance-covariance matrix;
S(j,k)=sum(t $ (sp(t,j) and sp(t,k)), (sp(t,j)-msp(j))*(sp(t,k)-msp(k)))/card(t);
display s;

Parameter b(n,i) resource availability;
b(n,"land")=sum(j, ha(n,j));
display b;

Parameter A(n,i,j);
A(n,"land",j)$(xbar(n,j) gt 0)=ha(n,j)/xbar(n,j);
display A;




VARIABLES u(n,j)        farm deviations
          y(n,i)        resource shadow prices
          lambda(n,j)   activity shadow prices
          alpha(n)      absolute risk aversion coefficient
          Pq(j,k)       factorisation matrix
          D(j,k)        Diagonal matrix of Cholesky decomposition
          L(j,k)        Lower triangular matrix of Cholesky decomposition
          LS            obj function value;

POSITIVE VARIABLES y, lambda, D
alpha
;

EQUATIONS FOB        objective function
          VEQ(n,j)   economic equilibrium constraint
          VSQ(n,j)   constraint relating to the Q matrix
          DCHOK(j,k) Cholesky decomposition
          VCT(n)     constraint on farm total cost
          V_land(n,i)


;



FOB..            LS =E=  0.001*(0.5*sum((n,j), sqr(u(n,j)))+ sum((n,i), y(n,i)*b(n,i))
                       + sum((n,j), c(n,j)*(xbar(n,j)))+sum((n,j), lambda(n,j)*(xbar(n,j)+epsilon(n)))
                       + sum(n, ((alpha(n))*sum((j,k), (xbar(n,j))*S(j,k)*(xbar(n,k)))))
                       - sum((n,j), p(n,j)*(xbar(n,j))));

VEQ(n,j)..            c(n,j)
                      +alpha(n)*sum(k, S(j,k)*(xbar(n,k)))
                      + sum(i,A(n,i,j)*y(n,i))+lambda(n,j)
                       =G= p(n,j);

VSQ(n,j)..            c(n,j)+lambda(n,j) =E=  sum(k, sum(jj, Pq(j,jj)*Pq(k,jj))* (xbar(n,k))) + u(n,j);

VCT(n)..              sum(j, c(n,j)*xbar(n,j))=L=0.5*(sum (j, xbar(n,j)*sum(k, sum(jj, Pq(j,jj)*Pq(k,jj))* (xbar(n,k)))))
                                                +sum(j, u(n,j)*xbar(n,j));

DCHOK(j,jj)..         Pq(j,jj) =E= sum(k, L(j,k)*((D(k,jj)+0.00000000000001)**0.5));

V_land(n,i)..        sum(j, A(n,i,j)*xbar(n,j))=L=1.0001*b(n,i);


D.FX(j,k)$(ord(j) ne ord(k))=0;
L.FX(j,k)$(ord(j) eq ord(k))=1;
L.FX(j,k)$(ord(j) lt ord(k))=0;


MODEL PRISK /ALL/;
veq.m(n,j)=xbar(n,j);
alpha.l(n)=0.00000001;
D.l(j,jj)$(ord(j) eq ord(jj))=0.001;
L.l(j,jj)$(ord(j) gt ord(jj))=0.00001;

execError = 0;

SOLVE PRISK using NLP minimizing LS;


parameter qmat(j,k)   square matrix of the non linear cost function
          prof(n)     farm profit
          ccal(n,j)   calibration error
          ac(n,j)     estimated average cost
          mc(n,j)     estimated marginal cost
          TC(n)       calculated farm Total costs
          TCc(n)      estimated  farm Total costs
          prof_ha(n)  profit per ha
          rev(n)      Total revenue
          ch_land(n)  check of land constraint
;

qmat(j,k)             = sum(jj, Pq.L(j,jj)*Pq.L(k,jj) ) ;

rev(n)               =sum(j, p(n,j)*xbar(n,j));

prof(n)               =sum(j, p(n,j)*xbar(n,j))
                      - 0.5*sum((j,k), xbar(n,j)*qmat(j,k)*xbar(n,k)) - sum(j, u.l(n,j)*xbar(n,j))
                     ;


ccal(n,j)               =veq.m(n,j)/xbar(n,j)*100-100;

ac(n,j)                 =0.5*sum(k, qmat(j,k)*xbar(n,k))+u.l(n,j);

mc(n,j)                 =sum(k, (qmat(j,k)*xbar(n,k)))+u.l(n,j);

TC(n)                   =sum(j, c(n,j)*xbar(n,j));

TCc(n)                  =0.5*(sum (j, xbar(n,j)*(sum (k, qmat(j,k)*xbar(n,k)))))+sum(j, u.l(n,j)*xbar(n,j));

prof_ha(n)              =prof(n)/sum(i,b(n,i));

ch_land(n)              = (sum(j, (veq.m(n,j)/r(n,j)))-sum(i,b(n,i)))/sum(i,b(n,i)) *100;


display xbar,veq.m,ccal, qmat,u.l,c,lambda.l,y.l,A,  prof, ac,rev,TC,TCc, alpha.l,  prof_ha, ch_land;

execError = 0;




***--------------  TEST CALIBRATED MODEL -------------------
Scalar rho inflation factor of the variance-covariance matrix of prices;
rho=1;


VARIABLES
x(n,j)    calibrated activity levels
ob        obj function value;

Positive variables
x;

EQUATIONS
FOBcalib         obj function
VSTRU(n,i)       resource constraint;

FOBcalib..             ob=E=sum(n, sum(j, p(n,j)*x(n,j)))
                            -0.5* sum (n, (sum(j, x(n,j)*(sum (k, qmat(j,k)*x(n,k))))))-sum(n, sum(j, u.l(n,j)*x(n,j)))
                            -0.5*rho*sum(n, (alpha.l(n)*sum((j,k), (x(n,j))*S(j,k)*(x(n,k))))) ;

VSTRU(n,i)..           sum(j, A(n,i,j)*x(n,j))=E=b(n,i);

x.l(n,j)=xbar(n,j);

execError = 0;

MODEL CALIB /FOBcalib, VSTRU/;

SOLVE CALIB using NLP maximising ob;

parameter ccal2(n,j)   calibration check
          ha0(n,j)     hectares allocated to each activity of the calibrated model
          share0(n)    farmland share under maize
          obF(n)       farm level obj function value
;

ccal2(n,j)          =x.l(n,j)/xbar(n,j)*100-100;

ha0(n,j)$x.l(n,j)   =x.l(n,j)/r(n,j);


share0(n)             =A(n,"land","Maize")*x.l(n,"Maize")/b(n,"land")*100;

obF(n)              =      sum(j, p(n,j)*x.l(n,j))
                            -0.5* sum(j, x.l(n,j)*(sum (k, qmat(j,k)*x.l(n,k))))-sum(j, u.l(n,j)*x.l(n,j))
                            -0.5*rho*(alpha.l(n)*sum((j,k), (x.l(n,j))*S(j,k)*(x.l(n,k)))) ;
display
x.l,ha, ha0, ccal2, vstru.m;

*execute_unload "results_cal"




***  ----------- SCENARIOS SIMULATIONS : CHANGE IN THE CROP PRICE VOLATILITY-------------------------
Set sc scenario /R1*R21/;

Parameter InflFactor(sc)  Variance-covariance matrix of prices in the simulation scenarios;
InflFactor(sc)=1+(ord(sc)-1)/20;
display Inflfactor;


parameter xS(n,j,sc)     results on activty level
          haS(n,j,sc)    results on hectares allocated to each activity
          change(n,j,sc) results on percentage change in hectare allocation to each activity compared to the baseline
          share(n,sc)    results on farmland share under maize
          shareb(n,j,sc) results on farmland share for each activity
          sum_xS(j,sc)
          sum_shareb(j,sc)
          Srho           saving InflFactor
          obTr(sc)       report obT
          obFr(n,sc)     report obF;

Srho=rho;

VARIABLES
xnew(n,j)      new activity level
obTnew         new overall objective function value
obFnew(n)      new farm level objective function value;

Positive variables
xnew;

EQUATIONS
FOBTnew           new overall objective function
FOBFnew(n)        new farm level objective function
VSTRUnew(n,i)
;

FOBTnew..    obTnew=E=           sum(n, sum(j, p(n,j)*xnew(n,j)))
                            -0.5*sum (n, (sum(j, xnew(n,j)*(sum (k, qmat(j,k)*xnew(n,k))))))-sum(n, sum(j, u.l(n,j)*xnew(n,j)))
                            -0.5*rho*sum(n, (alpha.l(n)*sum((j,k), (xnew(n,j))*S(j,k)*(xnew(n,k))))) ;

FOBFnew(n).. obFnew(n)=E=        sum(j, p(n,j)*xnew(n,j))
                            -0.5* (sum(j, xnew(n,j)*(sum (k, qmat(j,k)*xnew(n,k)))))-sum(j, u.l(n,j)*xnew(n,j))
                            -0.5*rho*(alpha.l(n)*sum((j,k), (xnew(n,j))*S(j,k)*(xnew(n,k)))) ; ;

VSTRUnew(n,i)..             sum(j, A(n,i,j)*xnew(n,j))=E=b(n,i);


xnew.l(n,j)=x.l(n,j);

MODEL SIMUL /FOBTnew, FOBFnew, VSTRUnew/;



LOOP (sc,
      rho=Srho;
      rho=InflFactor(sc);

execError = 0;

Solve SIMUL using NLP maximising obTnew;

xS(n,j,sc)                = xnew.l(n,j);
haS(n,j,sc)               = xnew.l(n,j)/r(n,j);
change(n,j,sc)$ha0(n,j)   = haS(n,j,sc)/ha0(n,j)*100-100;
share(n,sc)               = A(n,"land","Maize")*xnew.l(n,"Maize")/b(n,"land")*100;
shareb(n,j,sc)            = A(n,"land",j)*xnew.l(n,j)/b(n,"land")*100;
obTr(sc)                  = obTnew.l;
obFr(n,sc)                = obFnew.l(n);
);

sum_xS(j,sc) =   sum(n, xS(n,j,sc));
sum_shareb(j,sc)=sum(n, (xS(n,j,sc)*A(n,"land",j)))/(sum((n,i), b(n,i)))*100;

option xS:4:2:1
       haS:4:2:1
       change:4:2:1
       share:5:1:1
       shareb:5:2:1

display xs, haS,change, share, shareb, sum_shareb;

execute_unload "results_simul.gdx" ;




















