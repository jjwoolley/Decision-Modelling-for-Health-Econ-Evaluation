rm(list=ls())
library(heemod)
library(tidyverse)
s_names<-c('A',"B", 'C', 'D')

#two parts
#part A: "by hand" probability sensitivity analysis
#each run through gives one sensitivity analyse
#still need to add in combination therapy strategy
#part B: 100% heemod probability sensitivity analysis
#cannot follow the exact distributions from book
#not sure which shape parameters

#test



##### Assumptions #####
#transition probabilities
tpA2A<-1251/1734 #transition probability from A to A
tpA2B<-350/1734
tpA2C<-116/1734
tpA2D<-17/1734

tpB2B<-731/1258
tpB2C<-512/1258
tpB2D<-15/1258

tpC2C<-1312/1749
tpC2D<-437/1749

#costs (in british sterling)
dmca<-1701 #direct medical costs associated with state A
dmcb<-1774 
dmcc<-6948

ccca<-1055  #community care costs associated with state B
cccb<-1278
cccc<-2059

cAZT<-2278 #Zidovudine drug cost
cLam<-2087 #Lamivudine drug cost

#Other parameters
RR<- 0.509 #treatment effect of Combination treatment vs monotherapy 
#changes transition probabilities by RR

cDR<-0.06 #annual discount rate of costs
oDR <-0 #annual discount rate of benefits


#create matrix of all transitions from a state
#include deterministic mean, alpha, and beta
#using those, create random draw and then probabilistic transitions



#Probabilistic Sensitivity Analysis

#need to have random draws for all transitions FROM a single state
#then normalize each random draw by the sum of all random draws from that state
#this makes sure probabilities all add up to 1
#use a gamma distribution if alpha < 300 and normal distribution otherwise

tpA<-tibble(alpha=c(1251,350,116,17),beta=c(483,1384,1618,1717)) %>%
  mutate(dterm=alpha/(alpha+beta))%>%
  rowwise() %>%
  mutate(rd=if_else(alpha<300,rgamma(1,shape=alpha),rnorm(1,mean=alpha,sd=sqrt(beta))))

tpA_rd_sum<-sum(tpA$rd)

tpA_rd<-tpA %>%
  mutate(pb=rd/tpA_rd_sum)

#now do the same for transitions from State B

tpB<-tibble(alpha=c(731,512,15),beta=c(527,746,1243)) %>%
  mutate(dterm=alpha/(alpha+beta)) %>%
  rowwise() %>%
  mutate(rd=if_else(alpha<300,rgamma(1,shape=alpha,scale=1),rnorm(1,mean=alpha,sd=sqrt(beta)))) 

tpB_rd_sum<-sum(tpB$rd)

tpB_rd<-tpB %>%
  mutate(pb=rd/tpB_rd_sum)

#now for transitions from State C

tpC_rd<-tibble(alpha=c(1312,437),beta=c(437,1312)) %>%
  mutate(dterm=alpha/(alpha+beta)) %>%
  rowwise() %>%
  mutate(rd=rbeta(1,shape1=alpha,shape2=beta))%>%
  mutate(pb=rd)

#now make sure the transitions add up to 1 using method from the book
#tpC_rd$pb[1]=1-tpC_rd$pb[2]


###### COSTS ######
#direct medical costs and community care costs
#book gives assumption that standard error is equal to the mean of costs
#shape is alpha, scale is beta
dc_costs<-tibble(mean=c(1701,1774,6948,1055,1278,2059),
                 se=c(1701,1774,6948,1055,1278,2059)) %>%
  mutate(alpha=mean^2/se^2,
         beta=se) %>%
  rowwise()%>%
  mutate(pb=qgamma(runif(1),shape=alpha,scale=beta))



###### Other Parameters ######
RR_pb<-tibble(dterm=0.509) %>%
  mutate(ln_mean=log(dterm))%>%
  mutate(ln_se=(log(0.71)-log(0.365))/(1.96*2)) %>%
  mutate(pb=exp(rnorm(1,mean=ln_mean,sd=ln_se)))

ppmat_trans<-define_transition(
  tpA_rd[1,'pb'],tpA_rd[2,'pb'],tpA_rd[3,'pb'],tpA_rd[4,'pb'],
  0,tpB_rd[1,'pb'],tpB_rd[2,'pb'],tpB_rd[3,'pb'],
  0,0,C,tpC_rd[2,'pb'],
  0,0,0,1,
  state_names = s_names
)

s_A<-define_state(
  cost=dc_costs[1,'pb']+dc_costs[4,'pb']+cAZT,
  utility=1
)

s_B<-define_state(
  cost=dc_costs[2,'pb']+dc_costs[5,'pb']+cAZT,
  utility=1
)

s_C<-define_state(
  cost=dc_costs[3,'pb']+dc_costs[6,'pb']+cAZT,
  utility=1
)

s_D<-define_state(
  cost=0,
  utility=0
)


strat_monotherapy<-define_strategy(
  transition=ppmat_trans,
  A=s_A,
  B=s_B,
  C=s_C,
  D=s_D
)

res_mod<-run_model(
  method="beginning",
  strat_monotherapy,
  init=c(1,0,0,0),
  cycles=20,
  cost=cost,
  effect=utility
)

summary(res_mod)

get_counts(res_mod)
























######## HEEMOD PSA ############
rm(list=ls())
library(heemod)
library(tidyverse)

par_mod <- define_parameters(
  #need to adjust rr
  x=0.509,
  rr = ifelse(markov_cycle <= 2, x, 1),
  cost_lami = ifelse(markov_cycle <= 2, 2086.5, 0), #Lamivudine drug cost
  cost_zido = 2278, #Zidovudine drug cost
  
  ##### Assumptions #####
  #transition probabilities
  tpA2A=1251/1734, #transition probability from A to A
  tpA2B=350/1734,
  tpA2C=116/1734,
  tpA2D=17/1734,
  
  tpA_sum=tpA2A+tpA2B+tpA2C+tpA2D,
  
  tpA_rr_sum=tpA2A+(tpA2B+tpA2C+tpA2D)*rr,
  
  #these equations make sure each row adds up to one during sensitivity analysis
  tpAA=tpA2A/tpA_sum,
  tpAB=tpA2B/tpA_sum,
  tpAC=tpA2C/tpA_sum,
  tpAD=tpA2D/tpA_sum,
  
  tpB2B=731/1258,
  tpB2C=512/1258,
  tpB2D=15/1258,
  
  tpB_sum=tpB2B+tpB2C+tpB2D,
  
  tpBB=tpB2B/tpB_sum,
  tpBC=tpB2C/tpB_sum,
  tpBD=tpB2D/tpB_sum,
  
  tpC2C=1312/1749,
  tpC2D=437/1749,
  
  tpC_sum=tpC2C+tpC2D,
  
  tpCC=tpC2C/tpC_sum,
  tpCD=tpC2D/tpC_sum,
  
  #costs (in british sterling)
  dmca=1701, #direct medical costs associated with state A
  dmcb=1774, 
  dmcc=6948,
  
  ccca=1055,  #community care costs associated with state B
  cccb=1278,
  cccc=2059,
  
  cDR=0.06, #annual discount rate of costs
  oDR =0 #annual discount rate of benefits
)

mat_mono <- define_transition(
  tpAA, tpAB, tpAC, tpAD,
  0,    tpBB, tpBC, tpBD,
  0,    0,    tpCC, tpCD,
  0,    0,    0,    1
)


mat_comb <- define_transition(
  C,  tpAB*rr,  tpAC*rr,       tpAD*rr,
  0,        C,  tpBC*rr,       tpBD*rr,
  0,        0,        C,       tpCD*rr,
  0,        0,        0,       1
)


A_mono <- define_state(
  cost_health = dmca+ccca,
  cost_drugs = cost_zido,
  cost_total = discount(
    cost_health + cost_drugs, cDR, first = T),
  life_year = 1
)
B_mono <- define_state(
  cost_health = dmcb+cccb,
  cost_drugs = cost_zido,
  cost_total = discount(
    cost_health + cost_drugs, cDR, first = T),
  life_year = 1
)
C_mono <- define_state(
  cost_health = dmcc+cccc,
  cost_drugs = cost_zido,
  cost_total = discount(
    cost_health + cost_drugs, cDR, first = T),
  life_year = 1
)
D_mono <- define_state(
  cost_health = 0,
  cost_drugs = 0,
  cost_total = discount(
    cost_health + cost_drugs, cDR, first = T),
  life_year = 0
)

A_comb <- define_state(
  cost_health = dmca+ccca,
  cost_drugs = cost_zido + cost_lami,
  cost_total = discount(
    cost_health + cost_drugs, cDR, first = T),
  life_year = 1
)
B_comb <- define_state(
  cost_health = dmcb+cccb,
  cost_drugs = cost_zido + cost_lami,
  cost_total = discount(
    cost_health + cost_drugs, cDR, first = T),
  life_year = 1
)
C_comb <- define_state(
  cost_health = dmcc+cccc,
  cost_drugs = cost_zido + cost_lami,
  cost_total = discount(
    cost_health + cost_drugs, cDR, first = T),
  life_year = 1
)
D_comb <- define_state(
  cost_health = 0,
  cost_drugs = 0,
  cost_total = discount(
    cost_health + cost_drugs, cDR, first = T),
  life_year = 0
)

mod_mono <- define_strategy(
  transition = mat_mono,
  A_mono,
  B_mono,
  C_mono,
  D_mono
)
## No named state -> generating names.
mod_comb <- define_strategy(
  transition = mat_comb,
  A_comb,
  B_comb,
  C_comb,
  D_comb
)
## No named state -> generating names.
res_mod <- run_model(
  mono = mod_mono,
  comb = mod_comb,
  parameters = par_mod,
  cycles = 20,
  cost = cost_total,
  effect = life_year,
  method = "beginning",
  init = c(1, 0, 0, 0)
)


#net monetary benefit = threshold * effect_diff - cost_diff
summary(res_mod, threshold = c(50000, 100000, 150000))
summary(res_mod)

def_dsa<-define_dsa(
  tpA2B,0,.9,
  tpA2A,0,.9,
  dmca,100,10000,
  dmcb,100,10000,
  x,.25,.75
)

#need to change it so that it gives a random draw like in excel.


res_dsa<-run_dsa(res_mod,dsa=def_dsa)
summary(res_dsa)
plot(res_dsa, type = "difference", result='icer', strategy='comb', limits_by_bars = )

#should i use binomial or beta distribution
#what shape parameters to use?
#shape2 (sd) parameters were selected relatively arbitrarily
def_psa<-define_psa(
  
  tpA2A~binomial(0.721,20),
  tpA2B~binomial(0.202,20),
  tpA2C~binomial(0.067,20),
  tpA2D~binomial(0.01,20),
  
  tpB2B~binomial(0.581,20),
  tpB2C~binomial(0.407,20),
  tpB2D~binomial(0.012,20),
  
  tpC2C~binomial(0.75,20),
  tpC2D~binomial(0.25,20),
  
  dmca~gamma(1701,1701),
  dmcb~gamma(1774,1774),
  dmcc~gamma(6948,6948),
  
  ccca~gamma(1055,1055),
  cccb~gamma(1278,1278),
  cccc~gamma(2059,2059),
  
  rr~binomial(0.509,10)
)


res_psa<-run_psa(
  res_mod,psa=def_psa,N=1000
)

summary(res_psa)

#cost-efficiency plane
plot(res_psa,type = "ce", threshold = 5000)

#cost-effectiveness acceptability curves
plot(res_psa,type ="ac", max_wtp = 7000, log_scale = F)

#expected value of perfect information
plot(res_psa,type = "evpi", max_wtp = 7000, log_scale = F)

#covariance analysis on the results
plot(res_psa,type = "cov", diff = T, threshold = 50000)

plot(res_psa,type = "cov", diff = T, threshold = 100000)

plot(res_psa,type = "cov", diff = T, threshold = 150000)

xx<-res_psa$psa



