dt <- user(1)
initial(time) <- 0
update(time) <- (step+1) *dt
## Core equations for transitions between compartments:
update(S) <- S - n_SE1 - n_SE2 + n_VS + n_R1S + n_R2S + n_R12S - n_SV
update(V) <- V - n_VE1 - n_VE2 - n_VS + n_SV
update(E1) <- E1 + n_SE1 + n_VE1 - n_E1IA1
update(E2) <- E2 + n_SE2 + n_VE2 - n_E2IA2
# update(I1) <- I1 + n_E1I1 - n_I1D + n_I1R1
update(new_I1_cases) <- n_E1I1r + n_E1I1d
update(I1r) <- I1r + n_E1I1r - n_I1rR1
update(I1d) <- I1d + n_E1I1d - n_I1dD

update(A1) <- A1 + n_E1A1 - n_A1R1
#update(I2) <- I2 + n_E2I2 - n_I2D - n_I2R2
update(new_I2_cases) <- n_E2I2r + n_E2I2d
update(I2r) <- I2r + n_E2I2r - n_I2rR2
update(I2d) <- I2d + n_E2I2d - n_I2dD

update(A2) <- A2 + n_E2A2 - n_A2R2
update(R1) <- R1 + n_I1rR1 + n_A1R1 - n_R1E12 - n_R1S
update(R2) <- R2 + n_I2rR2 + n_A2R2 - n_R2E21 - n_R2S
update(E12) <- E12 + n_R1E12 - n_E12IA12
update(E21) <- E21 + n_R2E21 - n_E21IA21
#update(I12) <- I12 + n_E12I12 - n_I12D - n_I12R12

update(new_I12_cases) <- n_E12I12r + n_E12I12d
update(I12r) <- I12r + n_E12I12r - n_I12rR12
update(I12d) <- I12d + n_E12I12d - n_I12dD

update(A12) <- A12 + n_E12A12 - n_A12R12
#update(I21) <- I21 + n_E21I21 - n_I21D - n_I21R12

update(new_I21_cases) <- n_E21I21r + n_E21I21d
update(I21r) <- I21r + n_E21I21r - n_I21rR12
update(I21d) <- I21d + n_E21I21d - n_I21dD

update(A21) <- A21 + n_E21A21 - n_A21R12
update(R12) <- R12 + n_I12rR12 + n_A12R12 + n_I21rR12 + n_A21R12 - n_R12S
update(D) <- D + n_I1dD + n_I2dD + n_I12dD + n_I21dD

update(new_Deaths) <- n_I1dD + n_I2dD + n_I12dD + n_I21dD

p_SV <- 1 - exp(-(1 - epsilon_alpha) * rho)
n_SV <- rbinom(S, p_SV)
p_VS <- 1- exp(-alpha)
n_VS <- rbinom(V, p_VS)
p_SE1 <- 1 - exp(-(beta1 * I1 + beta2 *A1 + beta5 * I21 + beta6 * A21)/N)
n_SE1 <- rbinom(S, p_SE1)
p_SE2 <- 1- exp(-(beta3 * I2 + beta4 * A2 + beta7 * I12 + beta8 * A12)/N)
n_SE2 <- rbinom(S, p_SE2)
p_RS <- 1- exp(-eta)
n_R1S <- rbinom(R1, p_RS)
n_R2S <- rbinom(R2, p_RS)
n_R12S <- rbinom(R12, p_RS)

p_VE1 <- 1- exp(-((1 - epsilon_LV1) * beta1 * I1 + (1 - epsilon_LV1) * beta5 * I21 + (1 - epsilon_LVA1) * beta2 * A1 + (1 - epsilon_LVA1) * beta6 * A21)/N)
n_VE1 <- rbinom(V, p_VE1)
p_VE2 <- 1- exp(-((1 - epsilon_LV2) * beta3 * I2 + (1- epsilon_LV2) * beta7 * I12 + (1- epsilon_LVA2) * beta4 * A2 + (1- epsilon_LVA2) * beta8 * A12)/N)
n_VE2 <- rbinom(V, p_VE2)

p_E1IA1 <- 1 - exp(-omega)
n_E1IA1 <- rbinom(E1, p_E1IA1)

#n_E1IA1multi[] <- rmultinom(n_E1IA1, prob_symptomatic)
#prob_symptomatic[1] <- p
#prob_symptomatic[2] <- 1 - p
#dim(prob_symptomatic) <- 2
#dim(n_E1IA1multi) <- 2
#n_E1I1 <- n_E1IA1multi[1]
#n_E1A1 <- n_E1IA1multi[2]
n_E1I1 <- rbinom(n_E1IA1, p)
n_E1A1 <- rbinom(n_E1IA1 - n_E1I1, 1)

#n_I1rd[] <- rmultinom(n_E1I1, I1cfr)
#I1cfr[1] <- 1 - cfr1
#I1cfr[2] <- cfr1
#dim(I1cfr) <- 2
#dim(n_I1rd) <- 2
#n_E1I1r <- n_I1rd[1]
#n_E1I1d <- n_I1rd[2]
n_E1I1r <- rbinom(n_E1I1, 1 - cfr1)
n_E1I1d <- rbinom(n_E1I1 - n_E1I1r, 1)


p_E2IA2 <- 1 - exp(-omega)
n_E2IA2 <- rbinom(E2, p_E2IA2)

#n_E2IA2multi[] <- rmultinom(n_E2IA2, qprob_symptomatic)
#qprob_symptomatic[1] <- q
#qprob_symptomatic[2] <- 1-q
#dim(qprob_symptomatic) <- 2
#dim(n_E2IA2multi) <- 2
#n_E2I2 <- n_E2IA2multi[1]
#n_E2A2 <- n_E2IA2multi[2]
n_E2I2 <- rbinom(n_E2IA2, q)
n_E2A2 <- rbinom(n_E2IA2 - n_E2I2, 1)

#n_I2rd[] <- rmultinom(n_E2I2, I2cfr)
#I2cfr[1] <- 1 - cfr2
#I2cfr[2] <- cfr2
#dim(I2cfr) <- 2
#dim(n_I2rd) <- 2
#n_E2I2r <- n_I2rd[1]
#n_E2I2d <- n_I2rd[2]
n_E2I2r <- rbinom(n_E2I2, 1-cfr2)
n_E2I2d <- rbinom(n_E2I2 - n_E2I2r, 1)

p_I1dD <- 1- exp(-delta1)
n_I1dD <- rbinom(I1d, p_I1dD)
p_I1rR1 <- 1- exp(-gamma1)
n_I1rR1 <- rbinom(I1r, p_I1rR1)
p_A1R1 <- 1- exp(-gamma1)
n_A1R1 <- rbinom(A1, p_A1R1)
p_I2dD <- 1- exp(-delta2)
n_I2dD <- rbinom(I2d, p_I2dD)
p_I2rR2 <- 1- exp(-gamma2)
n_I2rR2 <- rbinom(I2r, p_I2rR2)
p_A2R2 <- 1- exp(-gamma2)
n_A2R2 <- rbinom(A2, p_A2R2)

p_R1E12 <- 1- exp(-((1 - epsilon_LR1) * beta3 * I2 + (1 - epsilon_LR1A) * beta4 * A2 + (1- epsilon_LR1) * beta7 * I12 + (1- epsilon_LR1A) * beta8 * A12)/N)
n_R1E12 <- rbinom(R1, p_R1E12)
p_R2E21 <- 1- exp(-((1 - epsilon_LR2) * beta1 * I1 + (1 - epsilon_LR2A) * beta2 * A1 + (1- epsilon_LR2) * beta5 * I21 + (1- epsilon_LR2A) * beta6 * A21)/N)
n_R2E21 <- rbinom(R2, p_R2E21)

p_E12IA12 <- 1 - exp(-omega)
n_E12IA12 <- rbinom(E12,p_E12IA12)

#n_E12IA12multi[] <- rmultinom(n_E12IA12, prob_reinfectionsymptomatic)
#prob_reinfectionsymptomatic[1] <- p
#prob_reinfectionsymptomatic[2] <- 1-p
#dim(prob_reinfectionsymptomatic) <- 2
#dim(n_E12IA12multi) <- 2
#n_E12I12 <- n_E12IA12multi[1]
#n_E12A12 <- n_E12IA12multi[2]
n_E12I12 <- rbinom(n_E12IA12, q)
n_E12A12 <- rbinom(n_E12IA12 - n_E12I12, 1)

#n_E12I12rd[] <- rmultinom(n_E12I12, I12cfr)
#I12cfr[1] <- 1 - cfr2
#I12cfr[2] <- cfr2
#dim(I12cfr) <- 2
#dim(n_E12I12rd) <- 2
#n_E12I12r <- n_E12I12rd[1]
#n_E12I12d <- n_E12I12rd[2]
n_E12I12r <- rbinom(n_E12I12, 1-cfr2)
n_E12I12d <- rbinom(n_E12I12 - n_E12I12r, 1)

p_E21IA21 <- 1 - exp(-omega)
n_E21IA21 <- rbinom(E21, p_E21IA21)

#n_E21IA21multi[] <- rmultinom(n_E21IA21, qprob_reinfectionsymptomatic)
#qprob_reinfectionsymptomatic[1] <- q
#qprob_reinfectionsymptomatic[2] <- 1-q
#dim(qprob_reinfectionsymptomatic) <- 2
#dim(n_E21IA21multi) <- 2
#n_E21I21 <- n_E21IA21multi[1]
#n_E21A21 <- n_E21IA21multi[2]
n_E21I21 <- rbinom(n_E21IA21, p)
n_E21A21 <- rbinom(n_E21IA21 - n_E21I21, 1)

#n_E21I21rd[] <- rmultinom(n_E21I21, I21cfr)
#I21cfr[1] <- 1- cfr1
#I21cfr[2] <- cfr1
#dim(I21cfr) <- 2
#dim(n_E21I21rd) <-2
#n_E21I21r <- n_E21I21rd[1]
#n_E21I21d <- n_E21I21rd[2]
n_E21I21r <- rbinom(n_E21I21, 1-cfr1)
n_E21I21d <- rbinom(n_E21I21 - n_E21I21r, 1)

p_I12dD <- 1 - exp(-delta12)
n_I12dD <- rbinom(I12d, p_I12dD)
p_I12rR12 <- 1- exp(-gamma12)
n_I12rR12 <- rbinom(I12r, p_I12rR12)
p_A12R12 <- 1- exp(-gamma12)
n_A12R12 <- rbinom(A12, p_A12R12)
p_I21dD <- 1- exp(-delta21)
n_I21dD <- rbinom(I21d, p_I21dD)
p_I21rR12 <- 1- exp(-gamma21)
n_I21rR12 <- rbinom(I21r, p_I21rR12)
p_A21R12 <- 1- exp(-gamma21)
n_A21R12 <- rbinom(A21, p_A21R12)

I1 <- I1r + I1d
I2 <- I2r + I2d
I12 <- I12r + I12d
I21 <- I21r + I21d

N <- S + V + E1 + E2 + I1 + A1 +I2 + A2 + R1 + R2 + E12 + E21 + I12 + A12 + I21 + A21 + R12 + D

initial(S) <- S_ini
initial(V) <- V_ini
initial(E1) <- E1_ini
initial(E2) <- E2_ini
#initial(I1) <- I1_ini
initial(new_I1_cases) <- 0
initial(I1r) <- I1r_ini
initial(I1d) <- I1d_ini

initial(A1) <- A1_ini
#initial(I2) <- I2_ini
initial(new_I2_cases) <- 0
initial(I2r) <- I2r_ini
initial(I2d) <- I2d_ini

initial(A2) <- A2_ini
initial(R1) <- R1_ini
initial(R2) <- R2_ini
initial(E12) <- E12_ini
initial(E21) <- E21_ini
#initial(I12) <- I12_ini
initial(new_I12_cases) <- 0
initial(I12r) <- I12r_ini
initial(I12d) <- I12d_ini
initial(A12) <- A12_ini
#initial(I21) <- I21_ini
initial(new_I21_cases) <- 0
initial(I21r) <- I21r_ini
initial(I21d) <- I21d_ini
initial(A21) <- A21_ini
initial(R12) <- R12_ini
initial(D) <- 0
initial(new_Deaths) <- 0

S_ini <- user(1059487)
V_ini <- user(6877200)
E1_ini <- user(93102)
E2_ini <- user(3879)
#I1_ini <- user(12987)
I1r_ini <- user(279742)
I1d_ini <- user(1972)
A1_ini <- user(0)
#I2_ini <- user(265)
I2r_ini <- user(11656)
I2d_ini <- user(82)
A2_ini <- user(0)
R1_ini <- user(0)
R2_ini <- user(0)
E12_ini <- user(116)
E21_ini <- user(5707)
#I12_ini <- user(303)
I12r_ini <- user(121)
I12d_ini <- user(1)
A12_ini <- user(0)
#I21_ini <- user(6)
I21r_ini <- user(2914)
I21d_ini <- user(21)
A21_ini <- user(0)
R12_ini <- user(0)

beta1 <- user(0.5)
beta2 <- user(0.08)
beta3 <- user(0.3)
beta4 <- user(0.0495)
beta5 <- user(0.5)
beta6 <- user(0.08)
beta7 <- user(0.3)
beta8 <- user(0.0495)
alpha <- user(0.0555)
eta <- user(0.0027)
epsilon_alpha <-  user(0.0862)
rho <- user(0)
epsilon_LV1 <- user(0.895)
epsilon_LVA1 <- user(0.895)
epsilon_LV2 <- user(0.519)
epsilon_LVA2 <- user(0.519)

omega <- user(0.25)
p <- user(0.15)
q <- user(0.15)
cfr1 <- user(0.0215)
cfr2 <- user(0.0113)
delta1 <- user(0.0018256)
gamma1 <- user(0.1)
delta2 <- user(0.00018256)
gamma2 <- user(0.1)
epsilon_LR1 <- user(0.519)
epsilon_LR1A <- user(0.519)
epsilon_LR2 <- user(0.519)
epsilon_LR2A <- user(0.519)
delta12 <- user(0.00018256)
gamma12 <- user(0.037)
delta21 <- user(0.0018256)
gamma21 <- user(0.03)