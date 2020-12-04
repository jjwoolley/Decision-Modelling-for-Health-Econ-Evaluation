#Decision Modelling for Health Economic Evaluation
#By Briggs, Blaxton, and Sculpher

rm(list=ls())


#Exercise 2.5

##### ASSUMPTIONS ######

#cycles take place annually
#utilities are all 1 except 0 in Death
#all individuals start in state A
#goes for 20 years from initial
n_t<-20


#states
n_states<-4
#A: 200 < cd4 < 500
#B: cd4 < 200
#C: AIDS
#D: Death
v_n<-c("A","B","C","D")

#treatments
#Monotherapy: Lamivudine in each period
#Combination: Zidovudine and Lamivudine for first two periods and then only Lamivudine



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
rr<- 0.509 #treatment effect of Combination treatment vs monotherapy 
    #changes transition probabilities by RR

cDR<-0.06 #annual discount rate of costs
oDR <-0 #annual discount rate of benefits

#state utilities vector
v_u<-c(A=1,B=1,C=1,D=0)

#state costs vector
v_c_Monotherapy<-c(A=dmca+ccca+cAZT,B=dmcb+cccb+cAZT,C=dmcc+cccc+cAZT,D=0)
v_c_Combination<-c(A=dmca+ccca+cAZT+cLam,B=dmcb+cccb+cAZT+cLam,C=dmcc+cccc+cAZT+cLam,D=0)




##### Matrix Setups #####

#initialize cohort trace matrix
v_s_init<-c(A=1,B=0,C=0,D=0)
mM<- matrix(0,
            nrow=(n_t+1), ncol=n_states,
            dimnames=list(0:n_t,v_n))
mM[1,]<-v_s_init




#initialize transition probability matrix
m_P<-matrix(0,
            nrow=n_states, ncol=n_states,
            dimnames=list(v_n,v_n))
#fill in matrix
m_P["A","A"]<-tpA2A
m_P["A","B"]<-tpA2B
m_P['A',"C"]<-tpA2C
m_P['A','D']<-tpA2D
m_P['B','B']<-tpB2B
m_P['B','C']<-tpB2C
m_P["B","D"]<-tpB2D
m_P['C','C']<-tpC2C
m_P['C','D']<-tpC2D
m_P['D','D']<-1


##### Solution #####
for(t in 1:n_t){
  mM[t + 1, ] <- mM[t, ] %*% m_P
}

mM_monotherapy<-mM[-1,]

#calculate costs
Costs_Monotherapy<-sum(mM[2:21,"A"])*(dmca+ccca+cAZT)+
  sum(mM[2:21,'B'])*(dmcb+cccb+cAZT) +
  sum(mM[2:21,'C'])*(dmcc+cccc+cAZT) 

yearly_costs_Monotherapy<-mM_monotherapy %*% v_c_Monotherapy
total_cost_Monotherapy<-sum(mM_monotherapy %*% v_c_Monotherapy)

#costs with discounting
v_dwc<-1/((1+cDR)^(1:(n_t)))

totalcost_discounted_Monotherapy<-sum(t(yearly_costs_Monotherapy) %*% v_dwc)

#life years
LY__Monotherapy<-sum(mM_monotherapy %*% v_u)



############# Combination ###########

#many values below taken as the same from monotherapy matrices




mM_comb<- matrix(0,
            nrow=(n_t+1), ncol=n_states,
            dimnames=list(0:n_t,v_n))
mM_comb[1,]<-v_s_init


m_P_comb_beg<-matrix(0,
            nrow=n_states, ncol=n_states,
            dimnames=list(v_n,v_n))

m_P_comb_beg["A","A"]<-1-tpA2B*rr-tpA2C*rr-tpA2D*rr
m_P_comb_beg["A","B"]<-tpA2B*rr
m_P_comb_beg['A',"C"]<-tpA2C*rr
m_P_comb_beg['A','D']<-tpA2D*rr
m_P_comb_beg['B','B']<-1-tpB2C*rr-tpB2D*rr
m_P_comb_beg['B','C']<-tpB2C*rr
m_P_comb_beg["B","D"]<-tpB2D*rr
m_P_comb_beg['C','C']<-1-tpC2D*rr
m_P_comb_beg['C','D']<-tpC2D*rr
m_P_comb_beg['D','D']<-1

for(t in 1:n_t){
  ifelse(t==1|t==2, mM_comb[t + 1, ] <- mM_comb[t, ] %*% m_P_comb_beg,
         mM_comb[t + 1, ] <- mM_comb[t, ] %*% m_P)
}

#delete initial state distribution before analysis
mM_Combination<-mM_comb[-1,]

#no discounting total cost
Costs_comb<-sum(mM_Combination[1:2,"A"])*(dmca+ccca+cAZT+cLam)+
  sum(mM_Combination[1:2,'B'])*(dmcb+cccb+cAZT+cLam) +
  sum(mM_Combination[1:2,'C'])*(dmcc+cccc+cAZT+cLam)+
  sum(mM_Combination[3:20,"A"])*(dmca+ccca+cAZT)+
  sum(mM_Combination[3:20,"B"])*(dmcb+cccb+cAZT)+
  sum(mM_Combination[3:20,"C"])*(dmcc+cccc+cAZT)

#yearly cost
yearly_costs_Monotherapy<-mM_monotherapy %*% v_c_Monotherapy
yearly_costs_comb<-cbind(c(mM_Combination[1:2,] %*% v_c_Combination,
                           mM_Combination[3:20,]%*%v_c_Monotherapy))

#costs with discounting
totalcost_discounted_comb<-sum(t(yearly_costs_comb) %*% v_dwc)

#life years
LY_Comb<-sum(mM_Combination %*% v_u)


#Incremental cost effectiveness ratio
ICER<-(totalcost_discounted_comb-totalcost_discounted_Monotherapy)/
  (LY_Comb-LY_Monotherapy)

################ HEEMOD ##################

library(heemod)
#copied from online
#https://cran.r-project.org/web/packages/heemod/vignettes/i_reproduction.html

par_mod <- define_parameters(
  rr = ifelse(markov_cycle <= 2, .509, 1),
  cost_lami = ifelse(markov_cycle <= 2, 2086.5, 0),
  cost_zido = 2278
)

mat_mono <- define_transition(
  1251/1734, 350/1734, 116/1734,  17/1734,
  0,         731/1258, 512/1258,  15/1258,
  0,         0,        1312/1749, 437/1749,
  0,         0,        0,         1.00
)
## No named state -> generating names.
mat_comb <- define_transition(
  C, 350/1734*rr, 116/1734*rr, 17/1734*rr,
  0, C,           512/1258*rr, 15/1258*rr,
  0, 0,           C,           437/1749*rr,
  0, 0,           0,           1.00)


A_mono <- define_state(
  cost_health = 2756,
  cost_drugs = cost_zido,
  cost_total = discount(
    cost_health + cost_drugs, .06, first = T),
  life_year = 1
)
B_mono <- define_state(
  cost_health = 3052,
  cost_drugs = cost_zido,
  cost_total = discount(
    cost_health + cost_drugs, .06, first = T),
  life_year = 1
)
C_mono <- define_state(
  cost_health = 9007,
  cost_drugs = cost_zido,
  cost_total = discount(
    cost_health + cost_drugs, .06, first = T),
  life_year = 1
)
D_mono <- define_state(
  cost_health = 0,
  cost_drugs = 0,
  cost_total = discount(
    cost_health + cost_drugs, .06, first = T),
  life_year = 0
)

A_comb <- define_state(
  cost_health = 2756,
  cost_drugs = cost_zido + cost_lami,
  cost_total = discount(
    cost_health + cost_drugs, .06, first = T),
  life_year = 1
)
B_comb <- define_state(
  cost_health = 3052,
  cost_drugs = cost_zido + cost_lami,
  cost_total = discount(
    cost_health + cost_drugs, .06, first = T),
  life_year = 1
)
C_comb <- define_state(
  cost_health = 9007,
  cost_drugs = cost_zido + cost_lami,
  cost_total = discount(
    cost_health + cost_drugs, .06, first = T),
  life_year = 1
)
D_comb <- define_state(
  cost_health = 0,
  cost_drugs = 0,
  cost_total = discount(
    cost_health + cost_drugs, .06, first = T),
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
summary(res_mod)


get_counts(res_mod,comb) 
