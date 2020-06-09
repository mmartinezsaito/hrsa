functions {
  // implicature encoding c: 1=simple, 2=complex, 11,21,33=unambiguous, 34=ambiguous
  // referent encoding r: 1=target, 2=competitor, 3=distractor 
  
  // literal listener
  vector UL0(int c) {
    vector[3] p;
    if (c == 1 || c == 2 || c == 34) p = [1, 1, 0]'/2.0;   
    else                             p = [1, 0, 0]'/1.0;  
    return p; 
  }  
  // exhaustifier
  vector UL1(int c) {
    vector[3] p;
    if      (c == 1)            p = [2, 1, 0]'/3.0;
    else if (c == 2 || c == 34) p = [1, 1, 0]'/2.0;   
    else                        p = [1, 0, 0]'/1.0;  
    return p; 
  }  
  // Gricean listener
  vector UL2(int c) {  
    vector[3] p;
    if      (c == 1)  p = [3, 1, 0]'/4.0;   
    else if (c == 2)  p = [3, 2, 0]'/5.0;   
    else if (c == 34) p = [1, 1, 0]'/2.0;   
    else              p = [1, 0, 0]'/1.0;   
    return p; 
  } 
}

data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1,upper=T> Ts[N];
  //int<lower=0,upper=4> trigger[N,T];
  //int<lower=1,upper=4> msg2[N,T];
  //int<lower=1,upper=4> msg3[N,T];
  //int<lower=1,upper=4> msg4[N,T];
  int<lower=0,upper=64> reft[N,T];
  int<lower=0,upper=64> refc[N,T];
  int<lower=0,upper=65> refd[N,T];
  int<lower=0,upper=3> choice[N,T]; 
  int<lower=0,upper=34> cond[N,T]; 
  matrix[82,8] sal;
}

transformed data {
  int nc = 3;  // number of conditions (two implicature types and unambiguous)
  int nr = 3;  // number of types of referents
  int nm = 4;  // number of messages
  int nch = 3; // number of referent tokens
  int nty = 3; // number of listener types
  matrix[82, 4] Sr = sal[:,1:4];
  matrix[82, 3] S = sal[:,6:8];
  row_vector[nty] typ = [0, 0, 1];   // listener type 2 
}

parameters {
  // to erase choice noise: delete 61 modify 122-3

  // group-level parameters
  real          mu_nc;
  real<lower=0> sd_nc;

  // fixed (yoked) parameter
  real<lower=0,upper=1> epsilon;        // uniform choice noise 
                                      // probability of choosing any referent no true of trigger message. subtracted from prob of the rest     

  // subject-level raw parameters, first follows norm(0,1), for later Matt Trick
  vector[N]    beta_nc;                
}

transformed parameters {
  
  // subject-level parameters
  vector<lower=0>[N] beta;         // inverse temperature

  // Matt Trick
  // the (group) hyperparameters Gaussian location and scale constrain the individual parameters
  for (i in 1:N) {
    beta[i] = exp( mu_nc + sd_nc * beta_nc[i] );  
  }

}

model {

  // ======= BAYESIAN PRIORS ======= //
  // group level hyperparameters
  mu_nc ~ normal(0, 1);
  sd_nc ~ cauchy(0, 1); // student_t(4,0,1); 
                        //cauchy(0,1);   // why half-Cauchy: Ahn, Haynes, Zhang (2017). 
                        //From v0.6.0, cauchy(0,5) -> cauchy(0,1) 

  // individual parameters: non-centered parameterization
  beta_nc ~ normal(0,1); // implies beta ~ exp(normal(mu_nc, sd_nc * beta_nc))

  // Franke, Degen (2016) fixing hyperparameters
  // beta ~ gamma(2, .5)   
  epsilon ~ gamma(.1, 2);                  // shape and rate parameters

 
  // ======= LIKELIHOOD FUNCTION ======= //
  // subject loop and trial loop

  for (i in 1:N) {

    for (t in 1:Ts[i]) {
      int c = cond[i, t];
      int ch = choice[i, t];
      vector[nch] u;
      vector[nch] p;
      int gridrow;
      real num0 = (c==1 || c==2 || c==34)? 1 : 2;
      real epsfac = num0 / (nch - num0);

      u =  typ[1] * UL0(c);
      u += typ[2] * UL1(c);
      for (j in 1:rows(Sr)) 
        if (Sr[j,1]==c && Sr[j,2]==reft[i,t] && Sr[j,3]==refc[i,t] && Sr[j,4]==refd[i,t]) 
          {gridrow = j; break;}
      u += typ[3] * UL2(c) .* S[gridrow]';
      
      p = softmax(beta[i] * u) + epsilon;
      p /= sum(p);
      ch ~ categorical(p); 
    }
  }
}

generated quantities {

  real llh[N, T];

  // For posterior predictive check
  real        choice_sim[N,T];
  vector[nch] p_pred[N, T];

  // set all posterior predictions to 0 (avoids NULL values)
  for (i in 1:N) for (t in 1:T) choice_sim[i, t] = -1;

  // subject and trial loops
  for (i in 1:N) {
    for (t in 1:Ts[i]) {
      int c = cond[i, t];
      int ch = choice[i, t];
      vector[nch] u;
      vector[nch] p;
      int gridrow;
      real num0 = (c==1 || c==2 || c==34)? 1 : 2;
      real epsfac = num0 / (nch - num0);

      u =  typ[1] * UL0(c);
      u += typ[2] * UL1(c);
      for (j in 1:rows(Sr)) 
        if (Sr[j,1]==c && Sr[j,2]==reft[i,t] && Sr[j,3]==refc[i,t] && Sr[j,4]==refd[i,t]) 
          {gridrow = j; break;}
      u += typ[3] * UL2(c) .* S[gridrow]';
      
      // ------- GQ ------- //
      p = softmax(beta[i] * u) + epsilon;
      p /= sum(p);
      llh[i, t] = categorical_lpmf(ch | p); 

      p_pred[i, t] = p; 
      choice_sim[i, t] = categorical_rng(p); // generate posterior prediction for current trial

      //p_pred = softmax(beta[i] * u);
      //for (r in 1:nch) {
      //  if (p[r]==0) p_pred[r] = epsilon;
      //  else         p_pred[r] -= epsilon*epsfac;
      //  if (p_pred[r]<0) p_pred[r] = 0;
      //}
      //p_pred /= sum(p_pred);
      //choice_sim[i, t] = categorical_rng(p_pred); // generate posterior prediction for current trial
      // ------- END GQ ------- //
    }
  }
}

