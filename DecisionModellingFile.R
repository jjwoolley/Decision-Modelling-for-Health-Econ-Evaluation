#Decision Modelling for Health Economic Evaluation
#By Briggs, Blaxton, and Sculpher




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
RR<- 0.509 #treatment effect of Combination treatment vs monotherapy 
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

#calculate life years 
LY_Monotherapy<-sum(mM[2:21,c('A','B','C')])

#calculate costs
Costs_Monotherapy<-sum(mM[2:21,"A"])*(dmca+ccca+cAZT)+
  sum(mM[2:21,'B'])*(dmcb+cccb+cAZT) +
  sum(mM[2:21,'C'])*(dmcc+cccc+cAZT) 

u<-sum(mM_monotherapy %*% v_u)
cc<-mM_monotherapy %*% v_c_Monotherapy
c<-sum(mM_monotherapy %*% v_c_Monotherapy)


#discounted costs
v_dwc<-1/((1+cDR)^(1:(n_t)))

sum(t(cc) %*% v_dwc)






################ HEEMOD ##################
library(heemod)
s_names<-c('A',"B", 'C', 'D')
mat_trans<-define_transition(
  tpA2A,tpA2B,tpA2C,tpA2D,
  0,tpB2B,tpB2C,tpB2D,
  0,0,tpC2C,tpC2D,
  0,0,0,1,
  state_names = s_names
)

s_A<-define_state(
  cost=dmca+ccca+cAZT,
  utility=1
)

s_B<-define_state(
  cost=dmcb+cccb+cAZT,
  utility=1
)

s_C<-define_state(
  cost=dmcc+cccc+cAZT,
  utility=1
)

s_D<-define_state(
  cost=0,
  utility=0
)


strat_monotherapy<-define_strategy(
  transition=mat_trans,
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
