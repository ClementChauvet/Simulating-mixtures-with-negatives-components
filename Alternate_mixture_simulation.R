#Ce programme simule le mélange d'un même nombre de lois normales positives et négatives 
#dans le cas ou il est possible de les appareiller 2 à 2 tel que a_i*f_i >=b_i*g_i


#Matrice de paramétres avec a chaque ligne une nouvelle loi. 
#La colonne 1 représente les moyennes, la 2 les ecarts types et la 3 les coefficients du mélange. Les n/2 premières lois sont celles de coefficients positif et les suivantes de coefficient négatif
x<- matrix(c(1,1,0.5,0,1,1.7,0,1,-0.7,1,1,-0.5), ncol=3, byrow=TRUE)
#Nombre de variable à simuler
n=1000

positive <- function(mu_1, sigma_1, alpha, mu_2, sigma_2, beta){
  
  if (sigma_1 != sigma_2){
    m<- -(sigma_2^2*mu_1 - sigma_1^2*mu_2)/(sigma_2^2-sigma_1^2)
    result<- (((exp((m-mu_1)^2/(2*sigma_1^2))-exp((m-mu_2)^2/(2*sigma_2^2)))<= (sigma_2*alpha)/(sigma_1*(-beta))) & (sigma_1 > sigma_2))
  }
  
  if (sigma_1 == sigma_2){
    result <- (mu_1 == mu_2)
    
    
  }
  
  return(result)
}

costMatrix <- function(mixture){
  
  n <- nrow(mixture)
  accept <- matrix(ncol=(n/2), nrow=(n/2))
  
  for (i in 1:(n/2)){
    
    for (j in ((n/2)+1):n){
      print(mixture)
      if(positive(mixture[i,1], mixture[i,2], mixture[i,3], mixture[j,1], mixture[j,2], mixture[j,3])==TRUE){
        
        accept[i,j-(n/2)] <- 1
      }
      else {
        
        accept[i,j-(n/2)] <- 0
      }
      
      
    }
  }
  
  
  return(accept) 
}

cost=costMatrix(x)
print(cost)



##Implementation de la Methode Hongroise par sergioreyblanco: https://github.com/sergioreyblanco/hungarian_method
solutionChecking <- function(matrix, solution, costMatrix){
  solutionAux <- solution
  foundSolutions <- list()
  flag_found = FALSE
  counter = 0
  previousDuplicate = -1
  duplicate = anyDuplicated(solutionAux[seq(2, 2*nrow(matrix), by = 2)])
  
  while(duplicate != 0){
    
    for(j in 1:ncol(matrix)){
      if(matrix[duplicate, j] == 0){
        counter = counter + 1
      }
      
      if(matrix[duplicate, j] == 0 && j != solutionAux[2*duplicate]){
        
        s<-solutionAux
        s[2*duplicate] = j
        flag_found = FALSE
        for(p in foundSolutions){
          paux <- as.vector(unlist(p))
          if(all(s == paux) == TRUE){
            flag_found = TRUE
            break
          }
        }
        if(flag_found == FALSE){
          solutionAux <- s
          
          foundSolutions[[length(foundSolutions)+1]] <- list(solutionAux)
          
          break
        }
      } 
      
    }
    
    if(anyDuplicated(solutionAux[seq(2, 2*nrow(matrix), by = 2)]) == 0){
      break
    }
    
    if(counter == 1){
      flag_found = TRUE
    }
    counter = 0
    
    if(flag_found == TRUE){
      
      for(i in 1:nrow(matrix)){
        for(j in 1:ncol(matrix)){
          
          if(matrix[i, j] == 0){
            s<-solutionAux
            s[2*(i-1)+1] <- i
            s[2*(i-1)+2] <- j
            flag_found = FALSE
            for(p in foundSolutions){
              paux<-as.vector(unlist(p))
              if(all(s == paux) == TRUE){
                flag_found = TRUE
                break
              }
            }
            if(flag_found == FALSE){
              solutionAux <- s
              break
            }
          }
        }
        if(flag_found == FALSE){
          break
        }
      }
      foundSolutions[[length(foundSolutions)+1]] <- list(solutionAux)
    }
    
    previousDuplicate = duplicate
    
    duplicate = anyDuplicated(solutionAux[seq(2, 2*nrow(matrix), by = 2)])
    
    if(previousDuplicate == duplicate){
      duplicate <- which(duplicated(solutionAux[seq(2, 2*nrow(matrix), by = 2)], fromLast = TRUE)==TRUE)[1]
    }
  }
  
  solutionValue = 0
  for(k in 1:nrow(matrix)){
    if(costMatrix[solutionAux[2*k-1],solutionAux[2*k]] != M){
      solutionValue = solutionValue + costMatrix[solutionAux[2*k-1],solutionAux[2*k]]
    }  
  }
  
  return(list("solution" = solutionAux, "cost" = solutionValue))
}

M <- 10000
Hungarian <- function(costMatrix){
  A <- costMatrix
  
  #Step 1: calculation fo reduced matrix
  
  #the minimum of each row is calculated and they are subtracted
  min = 10*M
  for(i in 1:nrow(A)){
    for(j in 1:ncol(A)){
      if(A[i, j] < min){
        min = A[i, j]
      }
    }
    A[i,] <- A[i,]-min
    
    min = 10*M
  }
  
  #the minimum of each column is calculated and they are subtracted
  min = 10*M
  for(j in 1:ncol(A)){
    for(i in 1:nrow(A)){
      if(A[i, j] < min){
        min = A[i, j]
      }
    }
    A[,j] <- A[,j]-min
    
    min = 10*M
  }
  
  while(TRUE){  
    #Step 2: zero strikethrough cross out  
    assignedRows <- list()
    assignedZeros <- list()
    alreadyCrossout <- list()
    for(i in 1:nrow(A)){
      for(j in 1:ncol(A)){
        if(A[i, j] == 0 && ((i-1)*ncol(A)+j) %in% alreadyCrossout == FALSE){
          assignedZeros <- c(assignedZeros, ((i-1)*ncol(A)+j))
          assignedRows <- c(assignedRows, i)
          if((i+1) < nrow(A)){
            for(k in (i+1):nrow(A)){
              if(A[k, j] == 0  && ((k-1)*ncol(A)+j) %in% alreadyCrossout == FALSE){
                alreadyCrossout<-c(alreadyCrossout, ((k-1)*ncol(A)+j))
              }
            }
          }
          
          break;
        }
      }
    }
    
    rowMarks <- list()
    columnMarks <- list()
    for(i in 1:nrow(A)){
      if(i %in% assignedRows == FALSE){
        rowMarks <- c(rowMarks, i)
      }
    }
    
    for(i in rowMarks){
      for(j in 1:ncol(A)){
        if(A[i, j] == 0 && j %in% columnMarks == FALSE){
          columnMarks <- c(columnMarks, j)
        }
      }
    }
    for(i in columnMarks){
      for(j in 1:nrow(A)){
        if(((j-1)*nrow(A)+i) %in% assignedZeros && j %in% rowMarks == FALSE){
          rowMarks <- c(rowMarks, j)
        }
      }
    }
    
    rowLines <- list()
    columnLines <- list()
    columnLines <- columnMarks
    for(i in 1:nrow(A)){
      if(i %in% rowMarks == FALSE){
        rowLines <- c(rowLines, i)
      }
    }
    
    if(length(rowLines) + length(columnLines) == nrow(A)){
      break;
    }else{
      #Step 3: reduced matrix is updated
      min = M
      for(i in 1:nrow(A)){
        for(j in 1:ncol(A)){
          if(A[i, j] < min && i %in% rowLines == FALSE && j %in% columnLines == FALSE){
            min = A[i, j]
          }
        }
      }
      for(i in 1:nrow(A)){
        for(j in 1:ncol(A)){
          if(i %in% rowLines == FALSE && j %in% columnLines == FALSE){
            A[i, j] <- A[i, j] - min
          } else if(i %in% rowLines == TRUE && j %in% columnLines == TRUE){
            A[i, j] <- A[i, j] + min
          } else if(i %in% rowLines == TRUE && j %in% columnLines == FALSE){
            A[i, j] <- A[i, j]
          } else if(i %in% rowLines == FALSE && j %in% columnLines == TRUE){
            A[i, j] <- A[i, j]
          }
        }
      }
    }
  }  
  
  #Step 4: selection of the solution
  solution <- c()
  for(i in 1:nrow(A)){
    for(j in 1:ncol(A)){
      if(A[i, j] == 0){
        solution <- c(solution, i, j)
        break
      }
    }
  }
  
  list <- solutionChecking(A, solution, costMatrix)
  
  return(list)
}
res=Hungarian(-cost)
assignment=c()
for(i in 1:(length(res$solution)/2)){
  assignment=append(assignment,res$solution[i*2]+nrow(x)/2)
}
print(assignment)
print(res)




rejection <- function( M, g,i,j,n,x) {
  naccepts <- 0
  result.sample <- rep(NA, n)
  
  while (naccepts < n) {
    y <- rnorm(1,x[i,1],x[i,2])
    u <- runif(1)
    
    if ( u <= x[i,3]*dnorm(y,x[i,1],x[i,2])+x[j,3]*dnorm(y,x[j,1],x[j,2]) / (M*dnorm(y,x[i,1],x[i,2])) ) {
      naccepts <- naccepts + 1
      result.sample[naccepts] = y
    }
  }
  
  result.sample
}

negative_mixture_sampling<-function(x,n, assignment){
  coeff=c(0)
  result=vector(,n)
  for(i in 1:(nrow(x)/2)){
    j=assignment[i]
    coeff=append(coeff,coeff[length(coeff)]+x[i,3]+x[j,3])#compute the cumulative sum of alpha-beta for each couple, it will be used as weights
  }
  U=runif(n)
  for(i in 1:(length(coeff)-1)){
    to_simulate=which(U >= coeff[i] & U<coeff[i+1])
    number_simulation=length(to_simulate)
    j=assignment[i]
    
    simulated=rejection(x[i,3],g,i,j,number_simulation,x)
    result[to_simulate]=simulated
  }
  
  return(result)
  
}
negative_mixture_sampling(x,n,assignment)