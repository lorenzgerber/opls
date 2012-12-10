##' Orthogonal Projection to Latent Structures (OPLS)
##'
##' \tabular{ll}{
##' Package: \tab opls\cr
##' Type: \tab Package\cr
##' Version: \tab 0.1\cr
##' Date: \tab 12-12-2011\cr
##' License: \tab GPL (>=2)\cr
##' LazyLoad: \tab yes \cr
##' }
##'
##' Projection to Latent Structures with Orthogonal siganl correction.
##' The package contains several additional functions for model
##' crossvalidation and evaluation.
##'
##' @name opls-package
##' @aliases opls-package
##' @docType package
##' @title Orthogonal Projection to Latent Structures (OPLS)
##' @author Lorenz Gerber \email{lorenz.gerber@@slu.se}
##' @keywords package
NA

##' Function to perform Orthogonal Projection to latent structures.
##' 
##' opls provides a simple function for calculating opls models as
##' described in Trygg and Wold (2001).
##' @export
##' @title opls
##' @param X N x K matrix, where N = observations and K = predictors
##' @param y response variable, dummy variable in case of DA class model 
##' @param nLVs number of latent variables to calculate, \code{numeric}
##' @param preprocess_method 0 = no preprocessing, 1 = centered, 2 = centered and 
##' @return object \code{opls_result}; a list containing:
##' \item{p_opls}{predictive X loadings}
##' \item{u_opls}{y score vector}
##' \item{w_opls}{predictive variable X weights}
##' \item{t_opls}{X score vector}
##' \item{c_opls}{y weights}
##' \item{b_opls}{y residuals}
##' \item{t_ortho}{orthogonal score vector(s)}
##' \item{p_ortho}{orthogonal X loadings}
##' \item{w_ortho}{orthogonal X weights}
##' \item{Eo_PLS}{X residuals}
##' \item{pre_vectors}{preprocessing X vector}
##' \item{avg_y}{preprocessing y vector}
##' @author Lorenz Gerber \email{lorenz.gerber@@slu.se}
opls<-function(X,y,nLVs,preprocess_method){

  ### Initiations

  T_ortho<-NULL    		### NULL vector
  P_ortho<-NULL			### NULL vector
  W_ortho<-NULL			### NULL vector

  X<-as.matrix(X)
  avg_y<-opls_preproc(t(as.matrix(y)),1)$vectors[1,]
  y<-t(t(y))

  if (preprocess_method==0){
    pre_data_vector<-opls_preproc(X,0)     ### If statements not finished!

  } else if(preprocess_method==1){
    pre_data_vector<-opls_preproc(X,1)

  } else if (preprocess_method==2){
    pre_data_vector<-opls_preproc(X,2)
  }

  avg_X<-pre_data_vector$vectors[1,]
  std_X<-pre_data_vector$vectors[2,]
  X<-pre_data_vector$data
  Eo_PLS<-X

  ##########################################################################
  ### Start algorithm
  #1
  w<-(t(y)%*%X)/as.vector((t(y)%*%y))		### Dot (scalar) product calculations
  w<-t(w)
  #2
  w<-w/as.vector(sqrt(t(w)%*%w))			# Generates vector
  ##########################################################################
  ### Start calculation for ortho factors first and predictive in the end
  #3
  for (counterLVs in 1:nLVs){
    te<-X%*%w/as.vector(t(w)%*%w)		# Generates vector
    #4
    ce<-(t(te)%*%y)/as.vector(t(te)%*%te)	# Generates scalar
    ce<-t(ce)
    #5
    u<-(y%*%ce)/as.vector(t(ce)%*%ce)		# Generates vector
    #6
    p<-t(te)%*%X/as.vector(t(te)%*%te)
    p<-t(p)
    ### End of calculation of predictive vector
    ### The algorithms will stop here after having already calculated the orthogonol components

    if (counterLVs<nLVs){
      #7
      w_ortho<-p-(as.vector(t(w)%*%p)/as.vector(t(w)%*%w))*(w)	
      #8
      w_ortho<-w_ortho/as.vector((sqrt(t(w_ortho)%*%w_ortho))) 
      #9
      t_ortho<-(X%*%w_ortho)/as.vector(t(w_ortho)%*%w_ortho) 
      #10
      p_ortho<-(t(t_ortho)%*%X)/as.vector(t(t_ortho)%*%t_ortho)	# Generates a row vector
      p_ortho<-t(p_ortho)
      #11
      Eo_PLS<-X-as.vector(t_ortho%*%t(p_ortho))
      #12
      T_ortho<-c(T_ortho,t_ortho)
      P_ortho<-c(P_ortho,p_ortho)
      W_ortho<-c(W_ortho,w_ortho)
      X<-Eo_PLS
    }
  }

  if(nLVs>1){
    T_ortho<-matrix(T_ortho,ncol=(nLVs-1)) ### reforming the vectors to matrices
  }

  if(nLVs>1){P_ortho<-matrix(P_ortho,ncol=(nLVs-1))}
  if(nLVs>1){W_ortho<-matrix(W_ortho,ncol=(nLVs-1))}

  #if(nLVs>1){T_ortho<-matrix(T_ortho,nrow=(nLVs-1),byrow=TRUE)}
  #if(nLVs>1){P_ortho<-matrix(P_ortho,nrow=(nLVs-1),byrow=TRUE)}
  #if(nLVs>1){W_ortho<-matrix(W_ortho,nrow=(nLVs-1),byrow=TRUE)}

  b<-w%*%ce
  
  #OPLS_result<-list(p_OPLS=p, u_OPLS=u,w_OPLS=w, t_OPLS=te, c_OPLS=ce, b_OPLS=b, T_ortho=T_ortho, P_ortho=P_ortho,W_ortho=W_ortho,Eo_PLS=X, pre_vectors=pre_data_vector$vectors, avg_y=avg_y)

  opls_result             <- list()
  opls_result$p_opls      <- p
  opls_result$u_opls      <- u
  opls_result$w_opls      <- w
  opls_result$t_opls      <- te
  opls_result$c_opls      <- ce
  opls_result$b_opls      <- b
  opls_result$t_ortho     <- T_ortho
  opls_result$p_ortho     <- P_ortho
  opls_result$w_ortho     <- W_ortho
  opls_result$Eo_PLS      <- X
  opls_result$pre_vectors <- pre_data_vector$vectors
  opls_result$avg_y       <- avg_y
  
  return(opls_result)
}


##' Predict new observations with an existing OPLS model
##'
##' Function to predict the response variable y by using an existing
##' \code{opls} model. The predictores K in X has to be of
##' the same length as the matrix X used for building the OPLS model.
##' @export
##' @title opls_pred
##' @param X_pred N x K matrix where N = observations to predict and K = predictors (same as
##' those used in \code{OPLS_result} from \code{opls}).
##' @param y unprocessed response vector Y for model set, also dummy variable in the case of DA classifications.
##' @param OPLS_result Return object from \code{OPLS} function.
##' @return object \code{OPLS_prediction_result}
##' \item{y_pred}{predicted response vector Y}
##' \item{t_pred}{predicted score vector X}
##' \item{t_pred}{predicted orthogonal scores of X}
##' \item{Eo_pls}{residuals of X}
##' @author Lorenz Gerber \email{lorenz.gerber@@slu.se}
opls_pred <-function(X_pred, y, opls_result){

  b           <- opls_result$b_opls
  w           <- opls_result$w_opls
  P_ortho     <- opls_result$p_ortho
  W_ortho     <- opls_result$w_ortho
  pre_vectors <- opls_result$pre_vectors
    
  ### Initiations
  T_ortho_pred<-NULL

  ### preprocessing the samples with OPLS preprocessing vectors
  if (mean(pre_vectors[1,])==0){
    X_filtered<-X_pred
  } else if (mean(pre_vectors[2,])==1){   # centering
    X_pred_centered<-X_pred-matrix(rbind(rep(pre_vectors[1,],times=nrow(X_pred))),ncol=ncol(X_pred),byrow=TRUE)
    X_pred<-X_pred_centered
    X_filtered<-X_pred
  } else if (mean(pre_vectors[2,])!=1){   # UV-scaling
    X_pred_centered<-X_pred-matrix(rbind(rep(pre_vectors[1,],times=nrow(X_pred))),ncol=ncol(X_pred),byrow=TRUE)
    X_pred_auto<-X_pred_centered/matrix(rbind(rep(pre_vectors[2,],times=nrow(X_pred))),ncol=ncol(X_pred),byrow=TRUE)
    X_pred<-X_pred_auto
    X_filtered<-X_pred 
  }

  #16
  if (!is.null(P_ortho)){
    X_filtered<-NULL

    for (line_nr in 1:nrow(X_pred)){
      t_new_ortho<-NULL
      x_new<-X_pred[line_nr,]

      for (i in 1:ncol(P_ortho)){
        t_new_ortho_x<-(x_new%*%W_ortho[,i])/((t(W_ortho[,i]))%*%W_ortho[,i]) # Line vector?
        #17
	t_new_ortho<-c(t_new_ortho, t_new_ortho_x)
	#18
	e_new_OPLS<-x_new-t_new_ortho_x*(P_ortho[,i])  ### NB - no dot product (multiplication with scalar)
        x_new<-e_new_OPLS
      }
      
      X_filtered<-rbind(X_filtered,x_new) 
      T_ortho_pred<-rbind(T_ortho_pred, t_new_ortho)
    }
  }

  t_pred<-X_filtered %*% w
  y_pred<-X_filtered %*% b+mean(y)

  opls_pred_result              <- list()
  opls_pred_result$y_pred       <- y_pred
  opls_pred_result$t_pred       <- t_pred
  opls_pred_result$t_ortho_pred <- T_ortho_pred
  opls_pred_result$Eo_pls       <- X_filtered

  return(opls_pred_result)                       
}


##' Preprocessing of Data for OPLS multivariate analysis
##'
##' This function is used internally from \code{opls} to preprocess both X and Y.
##' Data can either be centered or centered and unitvariate scaled.
##' @export
##' @title Preprocessing 
##' @param X Data matrix N x K, while N = observations and K = predictors
##' @param preprocess_method 0 = no preprocessing, 1 = centering, 2 = centering and unitvariate scaled
##' @return object \code{opls_preproc_result}; a list containing:
##' \item{data}{preprocessed data matrix}
##' \item{vectors}{matrix with column averages in row one and column standard deviations in row two.}
##' @author Lorenz Gerber \email{lorenz.gerber@@slu.se}
opls_preproc<-function (X,preprocess_method){
  if (preprocess_method==0){    ### No preprocessing
    preprocessing_data<-X
    preprocessing_vectors<-rbind(matrix(0,1,ncol(X)),matrix(1,1,ncol(X))) 

    opls_preproc_result         <- list()
    opls_preproc_result$data    <- preprocessing_data
    opls_preproc_result$vectors <- preprocessing_vectors
    
    return (opls_preproc_result)				### In R, a function can only return one varible

  } else if (preprocess_method==1){ 			### Centered (subtract mean)
    avg_X<-apply(X,2,mean) 			### Generates vector with mean for each column
    X_centered<-X-matrix(rbind(rep(avg_X,times=nrow(X))),ncol=ncol(X),byrow=TRUE)
    preprocessing_data<-X_centered
    preprocessing_vectors<-rbind(matrix(avg_X,1,ncol(X)),matrix(1,1,ncol(X))) 

    opls_preproc_result         <- list()
    opls_preproc_result$data    <- preprocessing_data
    opls_preproc_result$vectors <- preprocessing_vectors


    return (opls_preproc_result)

  } else if (preprocess_method==2){ 		### Centered and auto-scaled/UV-scaled (subtract mean and divied by SD)
    avg_X<-apply(X,2,mean)
    X_centered<-X-matrix(rep(avg_X,times=nrow(X)),ncol=ncol(X),byrow=TRUE)
    std_X<-apply(X_centered,2,sd)
    X_auto<-X_centered/matrix(rep(std_X,times=nrow(X)),ncol=ncol(X),byrow=TRUE)
    preprocessing_vectors<-rbind(avg_X,std_X)

    #temp<-as.double(sub(NaN,0,X_auto,fixed=TRUE))  ### if a Data variable has throughout zero values, this will result in NaN
    #X_auto<-matrix(temp,nrow=nrow(X_auto))	    ### therefore, here NaN is set back to zero, other all values in the OPLS function
    #rm(temp)					    ### becom NaN, one could actually also remove variables containing just zeroes

    preprocessing_data<-X_auto
    
    opls_preproc_result         <- list()
    opls_preproc_result$data    <- preprocessing_data
    opls_preproc_result$vectors <- preprocessing_vectors

    return (opls_preproc_result)
  }
}



##' Cross-validation for OPLS models 
##'
##' Cross validation by leave-one-out method. Yields R2 cumulative for X and Y
##' matrix. Further Q2 metrics predictive and orthogonal components are provided.
##' @export
##' @title opls_cv
##' @param X Matrix to be validated, N x K where N = observations and K = predictors
##' @param y Response vector.
##' @param nLVs Number of latent variables to be calculated.
##' @return object \code{output}; a list containing:
##' \item{PRESSp}{Predicted Sum of Squares of the predictive component}
##' \item{SSp}{Sum of Squares of the predictive component}
##' \item{Q2p}{Cross-validated R2 of the predictive component}
##' \item{R2xcum}{R2 cumulative of X}
##' \item{R2ycum}{R2 cumulative of y}
##' \item{PRESSo}{Predicted Sum of Squares of the orthogonal component(s)}
##' \item{SSo}{Sum of Squares of the orthogonal component(s)}
##' \item{Q2o}{Cross-validated R2 of the orthogonal component(s)}
##' @author Lorenz Gerber \email{lorenz.gerber@@slu.se}
opls_cv<-function(X,y,nLVs){
  leaveOutRow<-seq(1,nrow(X),1)

  fullmodel<-opls(X,y,nLVs,2)

  temp_PRESSp<-NULL
  temp_PRESSo<-NULL

  for(i in leaveOutRow){
    leaveOutOPLS<-opls(X[-i,],y[-i],nLVs,2)
    predictedleftout<-opls_pred(X,y,leaveOutOPLS)
    temp_PRESSp<-c(temp_PRESSp,(fullmodel$t_opls[i]-predictedleftout$t_pred[i])^2)
    if(nLVs>2)
      temp_PRESSo<-rbind(temp_PRESSo,(fullmodel$t_ortho[i,]-predictedleftout$t_ortho_pred[i,])^2)
    if(nLVs==1)
      temp_PRESSo<-c(temp_PRESSo,(fullmodel$t_ortho[i,]-predictedleftout$t_ortho_pred[i,])^2)
  }

  temp_PRESSp<-sum(temp_PRESSp)
  SSp<-sum(fullmodel$t_opls^2)
  Q2p<-1-temp_PRESSp/SSp

  if(nLVs>2)
    temp_PRESSo<-apply(temp_PRESSo,2,"sum")
  if(nLVs==2)
    temp_PRESSo<-sum(temp_PRESSo)

  if(nLVs>2)
    SSo<-apply(fullmodel$t_ortho^2,2,"sum")
  if(nLVs==2)
    SSo<-sum(fullmodel$t_ortho^2)
  if(nLVs>1)
    Q2o<-1-temp_PRESSo/SSo

  R2ycum<-1-sum(fullmodel$b_opls^2)/sum(((y-mean(y))/sd(y-(y-mean(y))))^2)
  R2xcum<-1-sum(apply(fullmodel$Eo_PLS^2,2,"sum"))/sum(apply((opls_preproc(X,2)$data)^2,2,"sum"))

  output          <- list()
  output$PRESSp   <- temp_PRESSp
  output$SSp      <- SSp
  output$Q2p      <- Q2p
  output$R2xcum   <- R2xcum
  output$R2ycum   <- R2ycum
  if(nLVs>1){
    output$PRESSo <- temp_PRESSo
    output$SSo    <- SSo
    output$Q2o    <- Q2o
  }
  

  return(output)

}



##' Quality metrics for a sequence of latent variables
##'
##' Function can be used to determine the ideal number of
##' components for an OPLS model.
##' @export
##' @title detcomponents
##' @param X Matrix to be evaluated, N x K where N = observations and K = predictors
##' @param y Response vector.
##' @param nLVseq Numbersequence for which quality metrics shall be calculated
##' @return Matrix with R2Xcum, R2ycum and Q2p 
##' @author Lorenz Gerber \email{lorenz.gerber@@slu.se}
detcomponents<-function(X,y,nLVseq){
  output<-matrix(rep(0,(3*length(nLVseq))),,3)


  for(i in nLVseq){
    OPLStemp<-OPLS_CV(X,y,i)
    output[i,]<-c(OPLStemp$R2xcum,OPLStemp$R2ycum,OPLStemp$Q2p)
  }

  colnames(output)<-c("R2Xcum","R2ycum","Q2")
  return(output)
}


##' Hotelling's T2 confidence calculation.
##' 
##' Calculates Hotelling´s T2 confidence borders for predictive
##' and first orthogonal component of an OPLS model.
##' @export
##' @title T2Hotelling
##' @param X Matrix to be evaluated, N x K where N = observations and K = predictors
##' @param y Response vector.
##' @param nLVs Number of latent variables to be calculated.
##' @param conf Confidence level.
##' @return object \code{output}; list containing:
##' \item{T2p}{Hotelling's T2 confidence border for the predictive component}
##' \item{T2o}{Hotelling's T2 confidence border for the first orthogonal component}
##' @author Lorenz Gerber \email{lorenz.gerber@@slu.se}
T2Hotelling<-function(X,y,nLVs,conf){
  tempOPLS<-OPLS(X,y,nLVs,2)
  T2p<-(apply(tempOPLS$t_opls,2,"sd")^2*qf(conf,nLVs,nrow(X)-nLVs)*nLVs*(nrow(X)^2-1)/(nrow(X)*(nrow(X)-nLVs)))^0.5
  T2o<-(apply(tempOPLS$t_ortho,2,"sd")[1]^2*qf(conf,nLVs,nrow(X)-nLVs)*nLVs*(nrow(X)^2-1)/(nrow(X)*(nrow(X)-nLVs)))^0.5

  output <- list()
  output$T2p <- T2p
  output$T2o <- T2o
  return(output)
}




##' Calculate ellipse points
##'
##' ellipse points,radially equispaced, given geometric par.s
##' @export
##' @title ellipsePoints
##' @param a length of half axes in x direction
##' @param b length of half axes in y direction
##' @param alpha angle (in degrees) for rotation
##' @param loc center of ellipse
##' @param n number of points
##' @return x,y coordinates of calculated ellipse
##' @author Martin Mächler
ellipsePoints<-function(a,b, alpha = 0, loc = c(0,0), n = 201){
  B <- min(a,b)
  A <- max(a,b)
  ## B <= A
  d2 <- (A-B)*(A+B)                   #= A^2 - B^2 
  phi <- 2*pi*seq(0,1, len = n)
  sp <- sin(phi)
  cp <- cos(phi)

  r <- a*b / sqrt(B^2 + d2 * sp^2)
  xy <- r * cbind(cp, sp)
  ## xy are the ellipse points for alpha = 0 and loc = (0,0)
  al <- alpha * pi/180
  ca <- cos(al)
  sa <- sin(al)
  xy %*% rbind(c(ca, sa), c(-sa, ca)) + cbind(rep(loc[1],n),rep(loc[2],n))
}
