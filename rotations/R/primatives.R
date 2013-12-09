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
  
  object<-matrix(object,dim(object))
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
  
  object<-matrix(object,dim(object))
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
  x<-matrix(x,dim(x))
  x<-x[i,...]
  return(as.Q4(x))
}

#' @S3method == Q4
#' @method == Q4
'==.Q4'<-function(e1,e2){

  e1<-matrix(e1,dim(e1))
  e2<-matrix(e2,dim(e2))
  if(all(e1==e2) || all(e1==-e2))  return(TRUE) else return(FALSE)

}

#' Arithmetic operators on SO(3)
#' 
#' These binary operators perform arithmetic on rotations in quaternion or rotation matrix form
#' (or objects which can be coerced into them).
#' 
#' The rotation group SO(3) is a multiplicative group so ``adding" rotations \eqn{R_1}{R1} and \eqn{R_2}{R2}
#' results in \eqn{R_1+R_2=R_2R_1}{R1+R2=R2R1}.  Similarly, the difference between rotations \eqn{R_1}{R1} and \eqn{R_2}{R2} is
#' \eqn{R_1-R_2=R_2^\top R_1}{R1-R2=R2'R1}.  With this definiton it is clear that 
#' \eqn{R_1+R_2-R_2=R_2^\top R_2R_1=R_1}{R1+R2-R2=R2'R2R1=R1}.  
#' If only one rotation is provided to subtraction then the inverse (transpose) it returned, 
#' e.g. \eqn{-R_2=R_2^\top}{-R2=R2'}.
#' 
#' @name Arithmetic
#' @aliases "+.SO3" "-.SO3" "+.Q4" "-.Q4"
#' @param x first arguement
#' @param y second arguement (optional for subtraction)
#' @return  \item{+}{the result of rotating the identity frame through x then y}
#'          \item{-}{the difference of the rotations, or the inverse rotation of only one arguement is provided}


NULL

#' @rdname Arithmetic
#' @aliases "-.SO3" "+.Q4" "-.Q4"
#' @S3method + SO3
#' @method + SO3

'+.SO3'<-function(x,y){

  y<-t(matrix(y,3,3))
  return(center(x,y))

}

#' @rdname Arithmetic
#' @aliases "+.SO3" "+.Q4" "-.Q4"
#' @S3method - SO3
#' @method - SO3

'-.SO3'<-function(x,y=NULL){

  if(is.null(y)) return(as.SO3(t(matrix(x,3,3))))
  
  return(center.SO3(x,y))
}

#' @rdname Arithmetic
#' @aliases "+.SO3" "-.SO3" "-.Q4"
#' @S3method + Q4
#' @method + Q4

'+.Q4'<-function(x,y){

  return(Q4(SO3(x)+SO3(y)))
}

#' @rdname Arithmetic
#' @aliases "+.SO3" "-.SO3" "+.Q4"
#' @S3method - Q4
#' @method - Q4

'-.Q4'<-function(x,y=NULL){
    
  if(is.null(y)){ 
    x[2:4]<- -1*x[2:4]
    return(x)
  }

  return(Q4(SO3(x)-SO3(y)))
}
