#' @S3method print SO3
#' @method print SO3

print.SO3<-function(x,...){
  Rs<-x
  len<-length(Rs)
  
  #if(len%%9!=0)
  #  stop("Input is not of the correct length.")
  
  if(len==9){
    print.default(as.SO3(matrix(Rs,3,3)),...)
  }else{
    print.default(Rs,...)
  }
}

#' @S3method head SO3
#' @method head SO3

head.SO3<-function(x,n=6L,...){
  
  #The following two lines are from 'head.matrix'
  stopifnot(length(n) == 1L)
  n <- if(n < 0L)  max(nrow(x) + n, 0L)  else min(n, nrow(x))
  
  x[seq_len(n) ,]
  
}

#' @S3method str SO3
#' @method str SO3

str.SO3<-function(object,...){
  
  n<-nrow(object)
  p<-ncol(object)
  object<-matrix(object,n,p)
  str(object)
}

#' @S3method print Q4
#' @method print Q4

print.Q4<-function(x,...){
  Qs<-x
  len<-length(Qs)
  
  #if(len%%4!=0)
  #  stop("Input is not of the correct length.")
  
  if(len==4){
    Qs<-matrix(Qs,1,4)
    negs<-length(which(Qs[2:4]<0))
    
    if(negs==0){ 
      
      print.default(bquote(.(Qs[1])+.(Qs[2])*i+.(Qs[3])*j+.(Qs[4])*k),...)
      
    }else if(negs==1){
      
      if(Qs[2]<0){
        Qs[2]<--Qs[2]
        print.default(bquote(.(Qs[1])-.(Qs[2])*i+.(Qs[3])*j+.(Qs[4])*k),...)
      }else if(Qs[3]<0){
        Qs[3]<--Qs[3]
        print.default(bquote(.(Qs[1])+.(Qs[2])*i-.(Qs[3])*j+.(Qs[4])*k),...)
      }else{
        Qs[4]<--Qs[4]
        print.default(bquote(.(Qs[1])+.(Qs[2])*i+.(Qs[3])*j-.(Qs[4])*k),...)
      }
      
    }else if(negs==2){
      
      if(all(Qs[2:3]<0)){
        
        Qs[2:3]<-abs(Qs[2:3])
        print.default(bquote(.(Qs[1])-.(Qs[2])*i-.(Qs[3])*j+.(Qs[4])*k),...)
        
      }else if(all(Qs[3:4]<0)){
        
        Qs[3:4]<-abs(Qs[3:4])
        print.default(bquote(.(Qs[1])+.(Qs[2])*i-.(Qs[3])*j-.(Qs[4])*k),...)
        
      }else{
        Qs[2:4]<-abs(Qs[2:4])
        print.default(bquote(.(Qs[1])-.(Qs[2])*i+.(Qs[3])*j-.(Qs[4])*k),...)
      }
    }else{ 
      Qs[2:4]<-abs(Qs[2:4])
      print.default(bquote(.(Qs[1])-.(Qs[2])*i-.(Qs[3])*j-.(Qs[4])*k),...)
    }
    
  }else{
    if(is.null(ncol(Qs))) {
      
      print.default(Qs,...)
      
    }else if(ncol(Qs)==4){
      
      colnames(Qs)<-c("Real","i","j","k")
      print.default(Qs,...)
      
    }else{
      
      print.default(Qs,...)
      
    }
  }
}

#print.Q4<-function(Qs,...){

#  len<-length(Qs)

#  if(len%%4!=0)
#    stop("Input is not of the correct length.")

#  if(len==4){
#    print.default(sprintf("%f + %f*i+ %f*j+%f*k ",Qs[1],Qs[2],Qs[3],Qs[4]),...)
#  }else{
#    colnames(Qs)<-c("Real","i","j","k")
#    print.default(Qs,...)
#  }
#}

#' @S3method head Q4
#' @method head Q4

head.Q4<-function(x,n=6L,...){
  
  #The following two lines are from 'head.matrix'
  stopifnot(length(n) == 1L)
  n <- if (n < 0L)  max(nrow(x) + n, 0L)  else min(n, nrow(x))
  
  x[seq_len(n), ]
  
}

#' @S3method str Q4
#' @method str Q4

str.Q4<-function(object,...){
  
  n<-nrow(object)
  p<-ncol(object)
  object<-matrix(object,n,p)
  str(object)
}

#' @S3method [ SO3
#' @method [ SO3
'[.SO3'<-function(x,i,...){
  x<-matrix(x,dim(x))
  x<-x[i,...]
  return(as.SO3(x))
}

#' @S3method [ Q4
#' @method [ Q4
'[.Q4'<-function(x,i,...){
  y<-matrix(x,length(x)/4,4)
  y<-y[i,...]
  return(as.Q4(y))
}

#' Arithmetic Operators on SO(3)
#' 
#' These binary operators perform arithmetic on rotations in quaternion or rotation matrix form
#' (or objects which can be coerced to them).
#' 
#' The rotation group SO(3) is a multiplicative group so addition of rotations \eqn{R1} and \eqn{R2}
#' is \eqn{R1+R2=R2R1}.  The difference between rotations \eqn{R1} and \eqn{R2} is
#' \eqn{R1-R2=R2'R1}.  With this definiton it is clear that \eqn{R1+R2-R2=R2'R2R1=R1}  Finally,
#' if only one rotation is provided to subtraction then the inverse (transpose) it returned, 
#' i.e. \eqn{-R2=R2'}.
#' 
#' @name Arithmetic
#' @aliases "+.SO3" "-.SO3" "+.Q4" "-.Q4"
#' @param x first arguement
#' @param y second arguement (optional for subtraction)
#' @return  \item{+}{the result of both rotations}
#'          \item{-}{the difference of the rotations, or the inverse rotation of only one arguement is provided}


NULL

#' @rdname Arithmetic
#' @aliases "-.SO3" "+.Q4" "-.Q4"
#' @S3method + SO3
#' @method + SO3

'+.SO3'<-function(x,y){
  x<-matrix(x,3,3)
  y<-matrix(y,3,3)
  return(as.SO3(y%*%x))
}

#' @rdname Arithmetic
#' @aliases "+.SO3" "+.Q4" "-.Q4"
#' @S3method - SO3
#' @method - SO3

'-.SO3'<-function(x,y=NULL){
  x<-matrix(x,3,3)
  if(is.null(y)) return(as.SO3(t(x)))
  y<-matrix(y,3,3)
  return(as.SO3(t(y)%*%x))
}

#' @rdname Arithmetic
#' @aliases "+.SO3" "-.SO3" "-.Q4"
#' @S3method + Q4
#' @method + Q4

'+.Q4'<-function(x,y){
  x<-SO3(x)
  y<-SO3(x)
  return(Q4(x+y))
}

#' @rdname Arithmetic
#' @aliases "+.SO3" "-.SO3" "+.Q4"
#' @S3method - Q4
#' @method - Q4

'-.Q4'<-function(x,y=NULL){
  
  if(is.null(y)){ 
    x[2:4]<--x[2:4]
    return(x)
  }
  x<-SO3(x)
  y<-SO3(y)
  return(Q4(x-y))
}
