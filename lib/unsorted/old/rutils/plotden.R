f
unction(x, y, z=NULL, cond=function(x) 1:length(x), title="", bw=1, rocfrac=1, nroc=2, adjust=0.5, ... ) {  
  if(!is.null(z)) {
    print('asldfkjhasldkfjlksd')                                        #r <- calc.roc(x,z,n=50)
    dx <- density(x,adjust=adjust)
    dy <- density(y,bw=dx$bw)    
    dz <- density(z,bw=dx$bw)
    ymax <- max(dx$y,dy$y,dz$y)
    plot(dx,ylim=c(0,ymax),main=paste(title),...)
    lines(dy,col=2)
    lines(dz,col=3)
  } else {
    r    <- getroc(x,y)
    dx   <- density(x)
    dy   <- density(y,bw=dx$bw)    
    ymax <- max(dx$y,dy$y)    
    xmx  <- quantile(c(dx$x,dy$x),0.5)
    plot(dx,ylim=c(0,ymax),xlim=c(0,xmx),main=paste(title,round(r,4)),...)
    lines(dy,col=2)    
    return(r)
  }  
}
