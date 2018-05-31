
Uniform_Prior = function(){
  
  mu = runif(1,min = 0, max = 1)
  return(mu)
}

Sine_Prior = function(){
  u = runif(1,min = 0 ,max = 1)
  mu = acos(1-2*u)/pi
  return(mu)
  
}

Cosine_Prior = function(){
  u = runif(1,min = 0 ,max = 1)
  f <- function(mu) {mu-1/pi* sin(pi*mu)-u}
  mu = uniroot(f,c(0,1))$root
  
}


CBT = function(n, prior, bn = log(log(n)), cn = log(log(n))){
  if(prior == "Uniform"){
    newmu <- Uniform_Prior
    mu_star = sqrt(2/n)
  }else if(prior == "Sine"){
    newmu <- Sine_Prior
    mu_star = (2*2*(2+1)/(pi^2*n))^(1/(2+1))
  }else if(prior == "Cosine"){
    newmu <- Cosine_Prior
    mu_star = (2*3*(3+1)/(pi^2*n))^(1/(3+1))
  }
      count = 0
      K = 1
      regret = 0 
      sum = 0 
      sum_sq = 0
      t = 0
      mu =  newmu()
      
      while(count < n){
        x = rgeom(1,prob = mu)
        if(count + x + 1 >= n){
          if(count + x + 1 == n){
            regret = regret + 1
            return(list("regret" = regret, "K" = K))
          }else{
            return(list("regret" = regret, "K" = K))
          }
        }
        
        count = count + x + 1
        regret = regret + 1
        t = t + x + 1
        sum = sum + 1
        sum_sq = sum_sq + 1
        
        x_bar = sum/t
        sigma = sqrt(sum_sq/t - (sum/t)^2)
        L = max(x_bar/bn , x_bar-cn*sigma/sqrt(t))
        if(L > mu_star){
          sum = 0 
          sum_sq = 0
          t = 0
          K = K + 1
          mu =  newmu()
        }
    }
  
}
  


Emp_CBT = function(n, prior, bn = log(log(n)), cn = log(log(n))){
  if(prior == "Uniform"){
    newmu <- Uniform_Prior
  }else if(prior == "Sine"){
    newmu <- Sine_Prior
  }else if(prior == "Cosine"){
    newmu <- Cosine_Prior
  }
  
  
    count = 0
    regret = 0 
    sum = 0 
    sum_sq = 0
    t = 0
    mu =  newmu()
    K = 1
   
    
      
        
      x = rgeom(1,prob = mu)
      if(x+1 >= n){
        if(x+1 == n){
          return(list("regret" = 1, "K" = 1))
        }else{
          return(list("regret" = 0, "K" = 1))
        }
      }
      
      count = count + x + 1
      regret = regret + 1
      t = t + x + 1
      sum = sum + 1
      sum_sq = sum_sq + 1
      
      indi = TRUE  
      wml = 1
      
      x_bar = sum/t
      sigma = sqrt(sum_sq/t - (sum/t)^2)
      L_min = max(x_bar/bn , x_bar-cn*sigma/sqrt(t))
      
      
      while(count < n){
        if(indi){
          x_bar = sum[K]/t[K]
          sigma = sqrt(sum_sq[K]/t[K] - (sum[K]/t[K])^2)
          L[K] = max(x_bar/bn , x_bar-cn*sigma/sqrt(t[K]))
          if(L[K]<L_min){
            wml = K
            L_min = L[K]
          }
        }else{
          x_bar = sum[wml]/t[wml]
          sigma = sqrt(sum_sq[wml]/t[wml] - (sum[wml]/t[wml])^2)
          L[wml] = max(x_bar/bn , x_bar-cn*sigma/sqrt(t[wml]))
          if(L[wml]<= L_min){
            L_min = L[wml]
          }else{
            wml = which.min(L)
            L_min = L[wml]
          }
        }
     
        
        if(L_min <= regret/n){
          indi = FALSE
          x = rgeom(1,prob = mu[wml])
          if(count + x + 1 >= n){
            if(count + x + 1 == n){
              regret = regret + 1
              return(list("regret" = regret, "K" = K))
            }else{
              return(list("regret" = regret, "K" = K))
            }
          }
            count = count + x + 1
            regret = regret + 1
            t[wml] = t[wml] + x + 1
            sum[wml] = sum[wml] + 1
            sum_sq[wml] = sum_sq[wml] + 1
          
        }else{
          
          indi = TRUE
          
          K = K +1
          mu[K] =  newmu()
          
          x = rgeom(1,prob = mu[K])
          if(count + x + 1 >= n){
            if(count + x + 1 == n){
              regret = regret + 1
              return(list("regret" = regret, "K" = K))
            }else{
              return(list("regret" = regret, "K" = K))
            }
          }
            count = count + x + 1
            regret = regret + 1
            t[K] =  x + 1
            sum[K] =  1
            sum_sq[K] =  1
          
        }
       
      }

}

Ana_CBT = function(n, data,  bn = log(log(n)), cn = log(log(n))){
  
  mK = ncol(data)
  
  ##shuffle the data
  data <- data[,sample(mK)]
  
  
  count = 0
  regret = 0 
  sum = 0 
  sum_sq = 0
  t = 0
  K = 1
  
  x = sample(data[,1],size = 1)
  
  count = count + 1
  regret = regret + x
  t = t + 1
  sum = sum + x
  sum_sq = sum_sq + x^2
  
  indi = TRUE  
  wml = 1
  
  x_bar = sum/t
  sigma = sqrt(sum_sq/t - (sum/t)^2)
  L_min = max(x_bar/bn , x_bar-cn*sigma/sqrt(t))
  
  
  while(count < n){
    if(indi){
      x_bar = sum[K]/t[K]
      sigma = sqrt(sum_sq[K]/t[K] - (sum[K]/t[K])^2)
      L[K] = max(x_bar/bn , x_bar-cn*sigma/sqrt(t[K]))
      if(L[K]<L_min){
        wml = K
        L_min = L[K]
      }
    }else{
      x_bar = sum[wml]/t[wml]
      sigma = sqrt(sum_sq[wml]/t[wml] - (sum[wml]/t[wml])^2)
      L[wml] = max(x_bar/bn , x_bar-cn*sigma/sqrt(t[wml]))
      if(L[wml]<= L_min){
        L_min = L[wml]
      }else{
        wml = which.min(L)
        L_min = L[wml]
      }
    }
    
    if(L_min <= regret/n){
      indi = FALSE
      x = sample(data[,wml],size = 1)

      count = count  + 1
      regret = regret + x
      t[wml] = t[wml] +  1
      sum[wml] = sum[wml] + x
      sum_sq[wml] = sum_sq[wml] + x^2
      
    }else{
      
      indi = TRUE
      
      K = K +1
      x = sample(data[,K],size = 1)
      
      count = count + 1
      regret = regret + x
      t[K] =    1
      sum[K] =  x
      sum_sq[K] =  x^2
      
    }
    
    if(K == mK){
      x_bar = sum[K]/t[K]
      sigma = sqrt(sum_sq[K]/t[K] - (sum[K]/t[K])^2)
      L[K] = max(x_bar/bn , x_bar-cn*sigma/sqrt(t[K]))
      if(L[K]<L_min){
        wml = K
        L_min = L[K]
      }
      
      for(i in (count+1):n){

          x = sample(data[,wml],size = 1)
         
          count = count  + 1
          regret = regret + x
          t[wml] = t[wml] +  1
          sum[wml] = sum[wml] + x
          sum_sq[wml] = sum_sq[wml] + x^2
        
          x_bar = sum[wml]/t[wml]
          sigma = sqrt(sum_sq[wml]/t[wml] - (sum[wml]/t[wml])^2)
          L[wml] = max(x_bar/bn , x_bar-cn*sigma/sqrt(t[wml]))
          if(L[wml]<= L_min){
            L_min = L[wml]
          }else{
            wml = which.min(L)
            L_min = L[wml]
          }
      }
      count = n
    }
      
      
  }
  
  return(list("regret" = regret, "K" = K))
}

