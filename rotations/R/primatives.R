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