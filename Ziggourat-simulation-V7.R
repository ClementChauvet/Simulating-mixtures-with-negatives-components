#Ce programme permet de simuler 100000 variables suivants la loi tronquée 2*N(2,2) -N(0,1)
rm(list=ls())
set.seed(5)

#Nombre de boites du ziggourat 
numberOfBoxes=30
#MAtrice définissant le mélange, la première colonne représente les coefficients, 
#la deuxième colonne représente les moyennes et le troisième représente les écarts types
X=matrix(c(2,2,2,-1,2,1),ncol=3,byrow=TRUE)
#Nombre de variables simulées
n=1000000
const=1/sqrt(2*pi)
f<-function(x,mu=0,sigma=1){
  return (const*exp((-(x-mu)^2)/(2*sigma^2))*(1/sigma))
}

ZiggouratInscrit=function(L,mu=0,sigma=1){
  #Genere un ziggourat sous la courbe
  SupMu=which(L>=mu)
  L1=L[SupMu]
  H=f(L1,mu,sigma)
  L2=L[-SupMu]
  H2=f(L2,mu,sigma)[1:length(L2)]
  print(L1)
  
  
  
  return(c(-H2,-H[2:length(H)]))
}
ZiggouratExte=function(L,mu=0,sigma=1){
  #Génere un ziggourat au dessus de la courbe
  SupMu=which(L>=mu)
  L1=L[SupMu]
  H=f(L1,mu,sigma)
  L2=c(L[-SupMu],mu)
  H2=f(L2,mu,sigma)[2:length(L2)]
  print(L1)
  return(c(H2,H[1:(length(H)-1)]))

  
}

tracer=function(L, hauteur,X){
  #Trace la courbe de densité du mélange X et le ziggourat définie par L et hauteur
  epsilon=0.00001
  Y=c()
  copie=c()
  Y[1]=L[1]
  for(i in 1:length(L)){
    Y[2*i]=L[i]
    Y[(2*i)+1]=L[i+1]-epsilon
    copie[2*i]=hauteur[i]
    copie[(2*i)+1]=hauteur[i]
    
  }
  plot(Y,copie,type="l")
  print((Y[length(Y)-1]-Y[2])/400)
  Y=seq(Y[2],Y[length(Y)-1], (Y[length(Y)-1]-Y[2])/400)
  Density=rep(0,length(Y))
  for(i in 1:nrow(X)){
    Density=Density+X[i,1]*dnorm(Y,X[i,2],X[i,3])
  }
  lines(Y,Density)
  
}


hauteur=list()


r1=qnorm(0.001,min(X[,2]),max(X[,3]))
r2=qnorm(0.999,max(X[,2]),max(X[,3]))
L=seq(r1,r2,abs(r1-r2)/numberOfBoxes)

L=c(X[,2],L)
L=unique(L)
L=sort(L)

#hauteur[[1]][i] représente la hauteur du ziggourat 1 entre points[[1]][i] et points[[1]][[i+1]]
for(i in 1:nrow(X)){
  if(X[i,1]>0){
    result=ZiggouratExte(L,X[i,2],X[i,3])
    hauteur<- c(hauteur,list(result*X[i,1]))
  }
  if(X[i,1]<0){
    result=ZiggouratInscrit(L,X[i,2],X[i,3])
    hauteur<- c(hauteur,list(result*-1*X[i,1]))
  }
  
}

H=rep(0,length(L)-1)
for(i in 1:length(hauteur)){
  H=H+hauteur[[i]]
}


tracer(L,H,X)



CalculAire=function(L,H){
  #Calcul l'aire sous le ziggourat
  Aire=0
  for(i in 1:(length(L)-1)){
    Aire=Aire+((L[i+1]-L[i])*H[i])
  }
  return(Aire)
}



CalculAireCumule=function(L,H){
  #Renvoie une liste de l'aire cumulé du ziggourat
  Aire=c(((L[2]-L[1])*H[1]))
  for(i in 2:(length(L)-1)){
    Aire=c(Aire,Aire[length(Aire)]+((L[i+1]-L[i])*H[i]))
  }
  
  return(Aire)
}



#simulate=function(L,H,X,r1,r2,n=1000){
  f=function(x){
    y=0
    for(i in 1:nrow(X)){
      y=y+X[i,1]*dnorm(x,X[i,2],X[i,3])
    }
    return(y)
  }
  
  #Calcul l'aire Cumule et la normalise
  AireCumule=CalculAireCumule(L,H)
  Aire=CalculAire(L,H)
  AireCumuleAjuste=AireCumule/Aire
  res=rep(0,n)
  nonGenere=1:n
  
  while(length(nonGenere)>0){
    U1=runif(length(nonGenere))
    nbAGenerer=rep(0,length(H))
    for(i in 1:length(H)){
      nbGenere=0
      #Choisis selon quel rectangle du ziggourat on va générer
      if(i==1){
        AGenerer=which(U1<AireCumuleAjuste[i])
      }
      else{
        AGenerer=which(U1<AireCumuleAjuste[i] & U1>=AireCumuleAjuste[i-1])
      }
      nbAGenerer[i]=length(AGenerer)
      
    
      #Une fois que l'on sait dans quel rectangle tirer Y, C'est une loi uniforme entre 2 points
      Y=runif((nbAGenerer[i]-nbGenere)*Aire,min=L[i],max=L[i+1])
      U2=runif((nbAGenerer[i]-nbGenere)*Aire)
      Acceptes=which(H[i]*U2<=f(Y))
      if(length(Acceptes)>0){
       res[nonGenere[AGenerer[1:min(length(Acceptes),nbAGenerer[i])]]]=Y[Acceptes[1:min(nbAGenerer[i],length(Acceptes))]]
        nbGenere=nbGenere+length(Acceptes)
      }
    
    }
    nonGenere=which(res==0)
  }
#  return(res)
#}
#res=simulate(L,H,X,r1,r2,n)
p1=hist(main="Répartition des valeurs générées",xlab="Valeurs",ylab="Répartition",res,breaks=100,freq=FALSE)
curve(2*dnorm(x,2,2)-dnorm(x,2,1),add=TRUE,col="red")





