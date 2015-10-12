eval_J<-function(P){
  #unpack parameters (design variables)
  a<-P[1:ndataset]
  b<-P[(ndataset+1):(2*ndataset)]
  alpha<-P[(2*ndataset+1)]
  beta<-P[(2*ndataset+2)]
  
  J = 0.0
  for(j in 1:ndataset){
    x<-as.vector(xydata[[j]][1])
    y<-as.vector(xydata[[j]][2])
    N<-length(x[,1])
    for(i in 1:N){
      sigma = alpha*abs(x[i,1])+beta
      J = J + log(sqrt(2*pi)*sigma)
      J = J + 0.5*((y[i,1]-a[j]*x[i,1]-b[j])^2)/(sigma^2)
    }
  }
  return(J)
}

eval_grad_J<-function(P){
  #unpack parameters (design variables)
  a<-P[1:ndataset]
  b<-P[(ndataset+1):(2*ndataset)]
  alpha<-P[(2*ndataset+1)]
  beta<-P[(2*ndataset+2)]
  
  dJdalpha = 0.0
  dJdbeta = 0.0
  dJda = c(rep(0.0,ndataset))
  dJdb = c(rep(0.0,ndataset))
  
  for(j in 1:ndataset){
    x<-as.vector(xydata[[j]][1])
    y<-as.vector(xydata[[j]][2])
    N<-length(x[,1])
    for(i in 1:N){
      #global parameters
      sigma = alpha*abs(x[i,1])+beta
      dsigmadalpha = abs(x[i,1])
      dJdalpha = dJdalpha + dsigmadalpha/(sigma) - abs(x[i,1])*((y[i,1]-a[j]*x[i,1]-b[j])^2)/(sigma^3)
      dJdbeta = dJdbeta + 1/(sigma) - ((y[i,1]-a[j]*x[i,1]-b[j])^2)/(sigma^3)
      #local parameters
      dJda[j] = dJda[j] - x[i,1]*(y[i,1]-a[j]*x[i,1]-b[j])/(sigma^2)
      dJdb[j] = dJdb[j] - (y[i,1]-a[j]*x[i,1]-b[j])/(sigma^2)
    }
  }
  
  dJ<-c()
  dJ[1:ndataset] = dJda
  dJ[(ndataset+1):(2*ndataset)] = dJdb
  dJ[(2*ndataset+1)] = dJdalpha
  dJ[(2*ndataset+2)] = dJdbeta
  
  return(dJ)
}

linear_estimator<-function(prefix,ndataset,P=c(rep(0,ndataset),rep(0,ndataset),1,1)){
  #make sure the optmization package is loaded
  library('nloptr')
  #read data into memory and store them as a list
  xydata<-list()
  fname<-c()
  for(j in 1:ndataset){
    fname[j]<-paste(prefix,j,".csv",sep="")
    xydata[[j]]<-read.csv(fname[j])
  }
  
  assign("xydata",xydata,envir=.GlobalEnv)
  assign("ndataset",ndataset,envir=.GlobalEnv)
  
  J<-eval_J(P)
  dJ<-eval_grad_J(P)
  
  res<-check.derivatives(.x=P,func=eval_J,func_grad=eval_grad_J,check_derivatives_print="all")
  
  res<-nloptr(x0=P,eval_f=eval_J,eval_grad_f=eval_grad_J,lb=c(rep(-Inf,(2*ndataset)),0,1e-10),ub=c(rep(+Inf,(2*ndataset)),+Inf,+Inf),opts=list("algorithm"="NLOPT_LD_LBFGS","print_level"=2,"check_derivatives"=TRUE,"xtol_rel"=1.0e-15,"maxeval"=500))
  
  print(res)
  #get solution
  P=as.vector(res[18])
  a<-P[[1]][1:ndataset]
  b<-P[[1]][(ndataset+1):(2*ndataset)]
  alpha<-P[[1]][(2*ndataset+1)]
  beta<-P[[1]][(2*ndataset+2)]
  
  #print optimality
  dJ<-eval_grad_J(P[[1]])
  J<-eval_J(P[[1]])
  normdJ<-norm(as.matrix(dJ))
  print(normdJ)
  print(dJ)
  print(J)
  
  #loop and plot data and estimator
  for(j in 1:ndataset){
    x<-as.vector(xydata[[j]][1])
    y<-as.vector(xydata[[j]][2])
    dev.new()
    #plot data
    plot(x[,1],y[,1])
    #plot linear model
    xplt=as.vector(range(x[,1]))
    yplt=a[j]*xplt+b[j]
    lines(xplt,yplt,type="l",col="red")
    #plot least squares
    lsfit<-lm(y[,1] ~ x[,1])
    f=as.vector(lsfit[1][1])
    bls=f[[1]][[1]]
    als=f[[1]][[2]]
    ypltls=als*xplt+bls
    lines(xplt,ypltls,type="l",col="blue")
    legend("topright",c("Model","Least-Squares"),lty=1,col=c("red","blue"))
    dev.copy(png,paste(prefix,j,".png",sep=""))
    dev.off()
  }
  
  data_out<-cbind(fname,a,b)
  write.csv(data_out,"solution.csv", row.names = FALSE,quote=FALSE)
  
  return(cbind(a,b,alpha,beta))
}


