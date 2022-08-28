rm(list=ls())

# Bibliotecas necessárias

library(quantmod)
library(Quandl)


# Carregamento dos dados 

getSymbols(c("^BVSP", "RENT3.SA","BRKM5.SA","SBSP3.SA"),
           periodicity='daily', 
           from='2018-12-01',
           to='2022-07-01'
)

dados <- merge(monthlyReturn(RENT3.SA[,6],type='log')[-1,], 
               monthlyReturn(BRKM5.SA[,6],type='log')[-1,], 
               monthlyReturn(SBSP3.SA[,6],type="log")[-1,])

names(dados) <- c( "RENT3","BRKM5","SBSP3")
Ibovespa = as.data.frame(merge(monthlyReturn(BVSP[,6],type='log')[-1,]))[,1]


m = 1 #      quantidade de restriçoes 
n = 3 # quantidade de icognitas de x 
a = c(1,1,1) # vetor com os coeficientes das restrições
b. = c(1) # vetor com os coeficientes do vetor independente
lambda = 1 # parâmetro de aversão ao risco (1/phi)
A = matrix(a, byrow = TRUE, nrow = m) # matriz dos coeficientes das restrições
b = matrix(b., byrow = FALSE, nrow = m) # vetor independente 
I = diag(1,n,n)
w = diag(1,m+n,1) # válido apenas para o caso de 1 restrição

# --------   Estimação do vetor de retorno esperado e
#             matriz de variâncias e covariâncias    ------ #


# método de estimação: 1- retornos históricos 
#                      2- modelo de índices 


Modelos <- function(ativo){
  return(summary(lm(ativo~Ibovespa)))
}

Pred <- function(ativo){
  return(lm(ativo~Ibovespa)$fitted.values)
}

Estimacao_Sigma <- function(Metodo_Estimacao){
  if (Metodo_Estimacao == 1){ 
    Sigma = cov(dados) 
    return(Sigma)
  }
  if(Metodo_Estimacao == 2){
    pred = c()
    for (i in 1:n){
      pred = as.data.frame(cbind(pred, Pred(dados[,i])))
      names(pred)[i] = names(dados)[i]
    }
    betas = c(rep(NA, n))
    for (i in 1:n){
      betas[i] <- Modelos(dados[,i])$coefficients[2]
    }
    Sigma <- matrix(c(rep(NA, n*n)), nrow = n) 
    for (i in 1:n) {
      for (j in 1:n){
        Sigma[i,j] <- betas[i]*betas[j]*var(Ibovespa)
      }
    }
    for (i in 1:n){
      Sigma[i,i] <- Sigma[i,i] + anova(lm(dados[,i]~Ibovespa))$`Mean Sq`[2]
    }
    return(Sigma)
  }
}

Estimacao_mu <- function(Metodo_Estimacao){
  if (Metodo_Estimacao == 1){
    return(apply(dados,2,mean))
  }
  else{
    pred = c()
    for (i in 1:n){
      pred = as.data.frame(cbind(pred, Pred(dados[,i])))
      names(pred)[i] = names(dados)[i]
    }
    return(apply(pred, 2, mean))
  }
}

Sigma = Estimacao_Sigma(2) # considerando modelo de índices
mu = Estimacao_mu(2) # considerando modelo de índices

# ---------  Método de Wolfe  --------- # 

q <-function(lambda, x, mu, Sigma){
  qx = lambda*t(x)*mu + (1/2)*t(x)*Sigma*x
  return(qx)
}


tableau <- function(A,Sigma,I,w,b,mu){
  frist = matrix(c(A,Sigma), byrow = TRUE, nrow = n+1)
  second = matrix(c(rep(0,n),-I),byrow = TRUE, nrow = n+1)
  third = matrix(c(rep(0,m),-t(A)),byrow = TRUE, nrow = n+1)
  fourth = matrix(c(rep(0,n),I), byrow = TRUE, nrow = n+1)
  fifth = matrix(c(rep(0,n),-I), byrow = TRUE, nrow = n+1)
  sixth = matrix(w,byrow = TRUE, nrow = n+1)
  finally = matrix(c(b,lambda*mu), byrow = TRUE, nrow = n+1)
  return(cbind(frist,second,third,fourth, fifth,sixth, finally))
}


# passo de inicialização do algoritmo


qx_inicial = c(rep(0,3*n+1),1)


if(sum(tableau(A,Sigma,I,w,b,mu)[,(4*n+3)] > 0) == m+n) {
  B = cbind(tableau(A, Sigma,I,w,b,mu)[,(2*n+2):(3*n+1)], tableau(A, Sigma,I,w,b,mu)[,4*n+2])
  NI = tableau(A, Sigma,I,w,b,mu)[,1:(2*n +1)]
  termo_independente = tableau(A, Sigma,I,w,b,mu)[,(4*n+3)]
  beta = c(rep(0,n),1)%*%solve(B)
  custo_relativo_inicial = c(qx_inicial[1:(2*n+1)] - beta%*%NI,rep(0,(n+1)))
  B <- rbind(B,c(rep(0,n), 1))
  NI<- rbind(NI, c(rep(0,2*n+1)))
  base <- B[-(m+n+1),]
  N = cbind(NI[,1:n], NI[,(2*n+1)])
  non_base <- N[-(m+n+1),]
  
}


if(sum(tableau(A,Sigma,I,w,b,mu)[,(4*n+3)] > 0) < m+n){
  B <- c()
  NI <- tableau(A, Sigma,I,w,b,mu)[,1:(2*n +1)]
  for(i in 1:n+m){
    if(tableau(A,Sigma,I,w,b,mu)[,(4*n+3)][i] < 0){
      B <- cbind(B,tableau(A, Sigma,I,w,b,mu)[,3*n+i])
    }
    if(tableau(A,Sigma,I,w,b,mu)[,(4*n+3)][i] >= 0) {
      B <- cbind(B,tableau(A, Sigma,I,w,b,mu)[,2*n+i])
    }
  }
  B <- cbind(B, w)
  termo_independente = tableau(A, Sigma,I,w,b,mu)[,(4*n+3)]
  beta = c(1,rep(0,n))%*%solve(B)
  custo_relativo_inicial = c(qx_inicial[1:(2*n+1)] - beta%*%NI,rep(0,n+1))
  B<- rbind(B,c(rep(0,n), 1))
  NI<- rbind(NI, c(rep(0,2*n+1)))
  N = cbind(NI[,1:n], NI[,2*n+1])
  base <- B[-(m+n+1),]
  non_base <- N[-(m+n+1),]
} 

basica = paste(seq(1,n+1),"*")
nao_basica = c(seq(1,2*n+1))
artificial = basica
originais = nao_basica

repeat{
  xb = solve(base)%*%termo_independente
  xn = c(rep(0,2*n+1))
  x = c(xn,xb)
  beta = B[(m+n+1),]%*%solve(base)
  custo_relativo_inicial = N[(m+n+1),] - beta%*%non_base 
  if(sum(custo_relativo_inicial < 0) == 0){
    cat("Solução ótima encontrada com y na base, problema original não é factível")
    break
  }
  col_pivo <- which.min(custo_relativo_inicial)
  auxiliar <- solve(base)%*%non_base[,col_pivo]
  if (sum(auxiliar > 0) == 0){
    cat("Problema é ilimitado") 
    break}
  ratio<-c()
  for(i in 1:(m+n)){
    if(auxiliar[i]<=0){
      ratio[i] <- NA}
    else{
      ratio[i] <- (xb[i]/auxiliar[i])}}
  lin_pivo = which.min(ratio) 
  x_out <- as.matrix(B[,lin_pivo])
  x_in <- as.matrix(N[,col_pivo])
  for (i in 1:dim(B)[2]){
    if(i == lin_pivo){
      B[,i] <- x_in
    }
  }
  for(i in 1:dim(N)[2]){
    if(i == col_pivo){
      N[,i] <- x_out
    }
  }
  ba = basica[lin_pivo]
  na = nao_basica[col_pivo]
  basica[lin_pivo] = na
  nao_basica[col_pivo] = ba
  base <- B[-(m+n+1),]
  non_base <- N[-(m+n+1),]
  if(sum(base[,dim(B)[2]] ==1)!= 1){
    cat("Soma dos wi's minimizada, pode seguir com a recursão")
    break
  }}

# recursao do método de wolfe

cb_z <- c()
for (i in 1:(n+1)){
  if (basica[i] %in% artificial){
    cb_z[i] <- 1}
  else{
    cb_z[i] <- 0}
}

rB <- rbind(base,cb_z)
rbase <- base
rN <- c()
rbasica <- basica
rnao_basica <- c()

non_base <- cbind(non_base[,1:n], tableau(A,Sigma,I,w,b,mu)[,(n+1):(2*n)],non_base[,ncol(non_base)])

for (i in 1:dim(non_base)[2]){
  if (nao_basica[i] %in% originais){
    rnao_basica[i] <- nao_basica[i]
    rN <- cbind(rN,non_base[,i])
  }
  rnao_basica <- c(na.omit(rnao_basica))}
rnon_base <- rN
rN <- rbind(rN,rep(0,dim(rN)[2]))


proximo <- function(rcol_pivo, rcustos){
  rcustos = rcustos[-rcol_pivo]
  pivo = which.min(rcustos)
  if(pivo < rcol_pivo){
    print(pivo)
  }
  else{
    print(pivo +1)
  }
}


repeat{
  if(sum(rbasica %in% artificial) > 0){
    rxb = solve(rbase)%*%termo_independente
    rxn = c(rep(0,dim(rN)[2]))
    rbeta = rB[(m+n+1),]%*%solve(rbase)
    rcusto_relativo_inicial = rN[(m+n+1),] - rbeta%*%rnon_base
    if(sum(rcusto_relativo_inicial < 0) == 0){
      cat("Solução ótima encontrada indices_básicas =", rbasica, "=> x =",solve(rbase)%*%termo_independente)
      break
    }  
    rcol_pivo <- which.min(rcusto_relativo_inicial) #verifica se é x ou v
    var = rnao_basica[rcol_pivo]
    if(var %in% originais[1:2*n]){
      if(var %in% originais[(n+1):(2*n)]){ # col_pivo é multipilcador de lagrange 
        corrs = which(var == originais)# testar se o x correspondente esta na base
        if (originais[corrs - n] %in% rbasica){ # v correspondente esta na base
          rcol_pivo = proximo(col_pivo, custos)
        }
        else{
          rcol_pivo = which.min(rcusto_relativo_inicial) 
        }}
      if(var %in% originais[(1):(n)]){ # col_pivo é variável x
        corrs = which(var == originais)# testar se o v correspondente esta na base
        if (originais[corrs + n] %in% rbasica){ # x correspondente esta na base
          rcol_pivo = proximo(rcol_pivo, rcusto_relativo_inicial)
        }
        else{
          rcol_pivo = which.min(rcusto_relativo_inicial) 
        }
        
      }
    }
    rauxiliar <- solve(rbase)%*%rnon_base[,rcol_pivo]
    if (sum(rauxiliar > 0) == 0){
      cat("Problema é ilimitado") 
      break}
    rratio<-c()
    for(i in 1:(m+n)){
      if(rauxiliar[i]<=0){
        rratio[i] <- NA}
      else{
        rratio[i] <- (rxb[i]/rauxiliar[i])}}
    rlin_pivo = which.min(rratio) 
    rx_out <- as.matrix(rB[,rlin_pivo])
    rx_in <- as.matrix(rN[,rcol_pivo])
    for (i in 1:dim(rB)[2]){
      if(i == rlin_pivo){
        rB[,i] <- rx_in
      }
    }
    for(i in 1:dim(rN)[2]){
      if(i == rcol_pivo){
        rN[,i] <- rx_out
      }
    }
    rba = rbasica[rlin_pivo]
    rna = rnao_basica[rcol_pivo]
    rbasica[rlin_pivo] = rna
    rnao_basica[rcol_pivo] = rba
    rbase <- rB[-(m+n+1),]
    rnon_base <- rN[-(m+n+1),]
    rxb = solve(rbase)%*%termo_independente
    if(sum(xb[which(rbasica %in% artificial)]) == 0){
      result = solve(rbase)%*%termo_independente 
      cat("Soma dos Zi's minimizada. Variáveis cujos valores são diferentes de zero: ",
          "indices(",rbasica[which(rbasica %in% seq(1,n))],") =", result[which(rbasica %in% seq(1,n))])
      break
    }
  }
}

