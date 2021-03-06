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
    if      (c == 1)            p = [1, 0, 0]';
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
  matrix[7,5] sal;
}

transformed data {
  int nc = 3;  // number of conditions (two implicature types and unambiguous)
  int nr = 3;  // number of types of referents
  int nm = 4;  // number of messages
  int nch = 3; // number of referent tokens
  int nty = 3; // number of listener types
  vector[nty] cpa0 = rep_vector(1, nty); // flat Dirichlet prior
  //matrix[nc, nr] S = [[1,1,1], [1,1,1], [1,1,1]] / 9.0; 
  vector[7]   Sr = sal[:, 1];
  matrix[7, 3] S = sal[:, 3:5];
}

parameters {
  // to erase choice noise: delete 61 modify 122-3


  // fixed (yoked) parameter
  //real<lower=0,upper=1> epsilon;        // uniform choice noise 
                                      // probability of choosing any referent no true of trigger message. subtracted from prob of the rest     
  real<lower=0> beta;         // inverse temperature

  // subject-level raw parameters, first follows norm(0,1), for later Matt Trick
  simplex[nty] typ[N];   // listener type (2 independent) probabilities
}

transformed parameters {


  // Matt Trick
  // the (group) hyperparameters Gaussian location and scale constrain the individual parameters

}

model {

  // ======= BAYESIAN PRIORS ======= //

  // individual parameters: centered parameterization
  // Franke, Degen (2016) fixing hyperparameters
  beta ~ gamma(2, .5);   
  //epsilon ~ gamma(.25, .1);                  // shape and rate parameters
  // typ ~ dirichlet(cpa0);                   // type distribution: fixed group hyperparameters
  for (i in 1:N) typ[i] ~ dirichlet(cpa0);     // type probabilities

 
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

      p  = typ[i, 1] * softmax(beta * UL0(c));
      p += typ[i, 2] * softmax(beta * UL1(c));
      for (j in 1:rows(Sr)) 
        if (Sr[j]==c) 
          {gridrow = j; break;}
      p += typ[i, 3] * UL2(c) .* S[gridrow]';
      
      //p += epsilon;
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
      vector[nch] p;
      int gridrow;
      real num0 = (c==1 || c==2 || c==34)? 1 : 2;
      real epsfac = num0 / (nch - num0);

      p  = typ[i, 1] * softmax(beta * UL0(c));
      p += typ[i, 2] * softmax(beta * UL1(c));
      for (j in 1:rows(Sr)) 
        if (Sr[j]==c) 
          {gridrow = j; break;}
      p += typ[i, 3] * UL2(c) .* S[gridrow]';

      // ------- GQ ------- //
      //p += epsilon;
      p /= sum(p);
      llh[i, t] = categorical_lpmf(ch | p); 

      p_pred[i, t] = p; 
      choice_sim[i, t] = categorical_rng(p); // generate posterior prediction for current trial
      // ------- END GQ ------- //
    }
  }
}
