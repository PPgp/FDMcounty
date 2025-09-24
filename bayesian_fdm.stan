// STAN model for estimating the FDM method
// Hana Sevcikova, 4/25/2025

// define the data objects
data {
  int<lower=0,upper=1> pre_working_age_in;        // 0 = no, 1 = yes
  int<lower=0,upper=1> working_age_in;            // 0 = no, 1 = yes
  int<lower=0,upper=1> retirement_in;             // 0 = no, 1 = yes
  int<lower=0,upper=1> post_retirement_in;        // 0 = no, 1 = yes
  int<lower=0,upper=1> pre_working_age_out;       // 0 = no, 1 = yes
  int<lower=0,upper=1> working_age_out;           // 0 = no, 1 = yes
  int<lower=0,upper=1> retirement_out;            // 0 = no, 1 = yes
  int<lower=0,upper=1> post_retirement_out;       // 0 = no, 1 = yes
  int<lower=0> N;                                 // number of age groups
  int<lower=0> Tm;                                // number of time periods
  vector[N] x;                                    // age vector
  array[Tm, N] int g;                             // age-specific net migrants (matrix of integers)
  matrix[Tm, N] pop;                              // age-specific population (matrix)
  matrix[Tm, N] pop_glob;                         // age-specific global population (matrix)
  vector[Tm] A;                                   // total in-migration counts for each year (vector)
  vector[Tm] B;                                   // total out-migration counts for each year (vector)
}

// Define parameters to be estimated, including restrictions on them
// Most of them are only defined if the corresponding RC component is included.
parameters {
  // in-RC parameters
  array[pre_working_age_in] real<lower=0, upper=1> alpha1_in;
  array[working_age_in] real<lower=0, upper=1> alpha2_in;
  array[retirement_in] real<lower=0, upper=1> alpha3_in;
  array[pre_working_age_in] real<lower=0, upper=1> a1_in;
  array[working_age_in] real<lower=a1_in, upper=1> a2_in;          // a2 >= a1
  array[retirement_in] real<lower=0, upper=1> a3_in;
  array[post_retirement_in] real<lower=0, upper=1> a4_in;
  array[working_age_in] real<lower=0,upper=55> mu2_in;             // 0 <= mu2 <= 55
  array[retirement_in] real<lower=55, upper=70> mu3_in;            // 55 <= mu3 <= 70
  array[working_age_in] real<lower=alpha2_in, upper=2> lambda2_in; // lambda2 >= alpha2
  array[retirement_in] real<lower=0, upper=2> lambda3_in;
  array[post_retirement_in] real<lower=-1,upper=0.05> lambda4_in;
  real<lower=0, upper=0.01> c_in;
  
  // out-RC parameters (the same restrictions as in-RC)
  array[pre_working_age_out] real<lower=0, upper=1> alpha1_out;
  array[working_age_out] real<lower=0, upper=1> alpha2_out;
  array[retirement_out] real<lower=0, upper=1> alpha3_out;
  array[pre_working_age_out] real<lower=0, upper=1> a1_out;
  array[working_age_out] real<lower=a1_out, upper=1> a2_out;
  array[retirement_out] real<lower=0, upper=1> a3_out;
  array[post_retirement_out] real<lower=0, upper=1> a4_out;
  array[working_age_out] real<lower=0, upper=55> mu2_out;
  array[retirement_out] real<lower=55, upper=70> mu3_out;
  array[working_age_out] real<lower=alpha2_out, upper=2> lambda2_out;
  array[retirement_out] real<lower=0, upper=2> lambda3_out;
  array[post_retirement_out] real<lower=-1,upper=0.05> lambda4_out;
  real<lower=0, upper=0.01> c_out;
  
  // the variance "v" parameter
  real<lower=0, upper=1> v;
}

transformed parameters {
  // derived parameters (visible to the outside)
  matrix[Tm, N] iota;     // matrix of time & age-specific in-migration counts
  matrix[Tm, N] omicron;  // matrix of time & age-specific out-migration counts
  matrix[Tm, N] delta;    // matrix of expected time & age-specific net-migration counts (i.e. the means)
  matrix[Tm, N] sigma;    // matrix of time & age-specific variances 
  vector[N] rx_in;        // age-specific RC curve for in-migration 
  vector[N] rx_out;       // age-specific RC curve for out-migration 
  matrix[Tm, N] rx_star_in;  // time & age-specific population-weighted RC curve for in-migration 
  matrix[Tm, N] rx_star_out; // time & age-specific population-weighted RC curve for out-migration 
  
  { // local variables (not visible outside of these parentheses)
    // they are used to compute the four individual RC components for in- and out- migration
    vector[N] mu_rc_1_in;
    vector[N] mu_rc_2_in;
    vector[N] mu_rc_3_in;
    vector[N] mu_rc_4_in;
    vector[N] mu_rc_1_out;
    vector[N] mu_rc_2_out;
    vector[N] mu_rc_3_out;
    vector[N] mu_rc_4_out;
    vector[N] zero;
    for(i in 1:N){
      zero[i] = 0;
    }
    // Compute the individual RC components and the totals;
    // if a component is missing, set it to a vector of 0 
    // in-migration
    mu_rc_1_in = pre_working_age_in==1?a1_in[1]*exp(-alpha1_in[1]*x):zero; 
    mu_rc_2_in = working_age_in==1?a2_in[1]*exp(-alpha2_in[1]*(x - mu2_in[1]) - exp(-lambda2_in[1]*(x - mu2_in[1]))):zero;
    mu_rc_3_in = retirement_in==1?a3_in[1]*exp(-alpha3_in[1]*(x - mu3_in[1]) - exp(-lambda3_in[1]*(x - mu3_in[1]))):zero;
    mu_rc_4_in = post_retirement_in==1?a4_in[1]*exp(lambda4_in[1]*(x)):zero;
    rx_in = mu_rc_1_in + mu_rc_2_in + mu_rc_3_in + mu_rc_4_in + c_in; // total RC curve for in-migration
    // out-migration
    mu_rc_1_out = pre_working_age_out==1?a1_out[1]*exp(-alpha1_out[1]*x):zero;
    mu_rc_2_out = working_age_out==1?a2_out[1]*exp(-alpha2_out[1]*(x - mu2_out[1]) - exp(-lambda2_out[1]*(x - mu2_out[1]))):zero;
    mu_rc_3_out = retirement_out==1?a3_out[1]*exp(-alpha3_out[1]*(x - mu3_out[1]) - exp(-lambda3_out[1]*(x - mu3_out[1]))):zero;
    mu_rc_4_out = post_retirement_out==1?a4_out[1]*exp(lambda4_out[1]*(x)):zero;
    rx_out = mu_rc_1_out + mu_rc_2_out + mu_rc_3_out + mu_rc_4_out + c_out; // total RC curve for out-migration
    
    // assure that the values are positive 
    for(i in 1:N){
      rx_in[i] = fmax(rx_in[i], 0);
      rx_out[i] = fmax(rx_out[i], 0);
    }
    
    for(j in 1:Tm){
      // compute population-weighted versions of the RC curves
      rx_star_in[j] = to_row_vector(rx_in) .* pop_glob[j]./sum(to_row_vector(rx_in) .* pop_glob[j]);
      rx_star_out[j] = to_row_vector(rx_out) .* pop[j]./sum(to_row_vector(rx_out) .* pop[j]);
      
      // compute age-specific in- and out-migration counts
      iota[j] = to_row_vector(rx_star_in[j] * A[j]); 
      omicron[j] = to_row_vector(rx_star_out[j] * B[j]);
      
      // compute the expected age-specfic net-migration
      delta[j] = iota[j] - omicron[j];
      
      // compute the variance
      sigma[j] = sqrt(iota[j]/v + omicron[j]/v);
    }
  } // end of the block that uses local variables
}


model {
  // likelihood
  for(j in 1:Tm){
    for(i in 1:N){
      g[j,i] ~ normal(delta[j,i], sigma[j,i]);
    }
  }
  //priors
  if(pre_working_age_in == 1){
    alpha1_in ~ normal(0,1);
    a1_in ~ normal(0,0.3);
  }
  if(pre_working_age_out == 1){
    alpha1_out ~ normal(0,1);
    a1_out ~ normal(0,0.3);
  }
  if(working_age_in == 1){
    alpha2_in ~ normal(0,1);
    a2_in ~ normal(0,0.3);
    mu2_in ~ normal(25,2);
    lambda2_in ~ normal(0,1);
  }
  if(working_age_out == 1){
    alpha2_out ~ normal(0,1);
    a2_out ~ normal(0,0.3);
    mu2_out ~ normal(25,2);
    lambda2_out ~ normal(0,1);
  }
  if(retirement_in == 1){
    alpha3_in ~ normal(0,1);
    a3_in ~ normal(0,0.3);
    mu3_in ~ normal(63,2);
    lambda3_in ~ normal(0,1);
  }
  if(retirement_out == 1){
    alpha3_out ~ normal(0,1);
    a3_out ~ normal(0,0.3);
    mu3_out ~ normal(63,2);
    lambda3_out ~ normal(0,1);
  }
  if(post_retirement_out == 1){
    a4_in ~ normal(0,0.05);
    lambda4_in ~ normal(0,0.01);
  }
  if(post_retirement_out == 1){
    a4_out ~ normal(0,0.05);
    lambda4_out ~ normal(0,0.01);
  }
  c_in ~ normal(0, 0.005);
  c_out ~ normal(0, 0.005);
  v ~ uniform(0, 1);
}


generated quantities{
  array[Tm, N] real g_pred;  // predictive distribution of net migration counts g
  {
    for(j in 1:Tm){
      for(i in 1:N){
        g_pred[j,i] = normal_rng(delta[j, i], sigma[j,i]);
      }
    }
  }
}
