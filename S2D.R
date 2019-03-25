# Source-to-dose module (S2D) for the Human Exposure Mofdel (HEM) 
# Written for EPA by Graham Glen at ICF, Feb - May 2017
# Last modified by GG on November 28, 2018
# This version adds chemical infiltration from outdoor air 
# Now all houses are run, even ones that use no products containing the chemical agent
wd <- "C:/main/HEM/Nov2018"
setwd(wd)

s2d = function(control.file="control_file.txt", number.of.houses= NULL) {
  library(data.table)
  library(stringr)
  library(plyr)
  library(dplyr)
  library(dtplyr)
  library(ggplot2)
  library(bit64)
  library(foreach)
  library(doParallel)
  library(reshape2)
  
  ###################################################################
  ############# Start of function definitions #######################
  ###################################################################
  
  # Distrib returns samples from distributions - Written by WGG for SHEDS-HT in 2012
  
  distrib = function(shape="",par1=NA,par2=NA,par3=NA,par4=NA,lt=NA,
                     ut=NA, resamp="y",n=1,q=NA,p=c(1),v="" ) {
    # Distrib generates samples from a distribution. If a vector of quantiles
    # q is given, the corresponding values are returned. Otherwise n determines 
    # the number of samples, but then quantiles are randomly generated first.
    # Written by Graham Glen, July 2012.
    
    m <- ""
    n <- round(n)
    if (n<=0) m <- paste("# samples requested = ",n)
    if (length(par1)>1|length(par2)>1|length(par3)>1|length(par4)>1) {
      m <- "Non-scalar distribution parameters"
    }
    if (is.na(q[1])) q <- runif(n)
    s <- strtrim(tolower(shape),4)
    r <- strtrim(tolower(resamp),1)
    if (is.na(lt)&is.na(ut)) r <- "n"
    if (!is.na(ut)&!is.na(lt)&ut<lt) m <- paste("Truncation limits",lt,"and",ut)
    if (!is.numeric(q)) m <- "Quantiles are not numeric"
    if (min(q)<0 || max(q)>1) m <- "Quantiles not all between 0 and 1" 
    
    if (m!="") {
    } else if (s=="bern" || s=="bino") {
      if (is.na(par1)) par1 <- 0.5
      lt <- NA 
      ut <- NA
      if (par1<0 | par1>1) { m <- paste("Invalid binomial parameter",par1)
      } else {
        x <- q
        x[q<=1-par1] <- 0
        x[q> 1-par1] <- 1
      }
    } else if (s=="beta") {
      if (is.na(par1)) par1 <- 1
      if (is.na(par2)) par2 <- 1
      if (is.na(par3)&&is.na(par4)) { par3 <- 0; par4 <- 1
      } else if (is.na(par4)) { par4 <- par3+1
      } else if (is.na(par3)) par3 <- par4-1
      if (!is.na(lt) && lt>par4) m <- "Lower beta truncation above par4"
      if (!is.na(ut) && ut<par3) m <- "Upper beta truncation below par3"
      if (!is.na(lt)&&!is.na(ut)&&lt>ut) m <- "Lower truncation above upper"
      w <- par4-par3
      if (par1<=0 | par2<=0 | w<0) { 
        m <- paste("Invalid beta parameters",par1,par2,par3,par4)
      } else if (r=="y") {
        if (!is.na(lt)) qlo <- pbeta((lt-par3)/w,par1,par2) else qlo <- 0
        if (!is.na(ut)) qhi <- pbeta((ut-par3)/w,par1,par2) else qhi <- 1
        q <- qlo+q*(qhi-qlo)
      } 
      if (m=="") x <- par3 + w*qbeta(q,par1,par2)
    } else if (s=="disc" || s=="prob") {
      if (mode(p)=="character") p <- scan(text=p,quiet=TRUE)
      if (mode(v)!="numeric") {
        if (mode(try(scan(text=v,quiet=TRUE),silent=TRUE))=="numeric") {
          v <- c(scan(text=v,quiet=TRUE),recursive=TRUE)
        } else {
          v <- c(scan(what=list(""),text=v,quiet=TRUE),recursive=TRUE)
        }  
      } 
      if (length(v)==0) v <- 0
      if (s=="prob") v <- 1:length(p)
      if (s=="disc" && length(p)==1 && length(v)>1) p <- rep(1,length(v))
      if (s=="disc" && length(v)==1 && length(p)>1) v <- 1:length(p)
      if (!is.na(lt)) lt <- min(v[which(v>=lt)])
      if (!is.na(ut)) ut <- max(v[which(v<=ut)])
      if (length(p)!=length(v)) { 
        m <- paste("Unequal vectors, v=",length(v),", p=",length(p))
      } else if (min(p)<0) { m <- "Negative probabilities found"  
      } else if (r=="y") {
        if (!is.na(lt)) p[v<lt] <- 0
        if (!is.na(ut)) p[v>ut] <- 0
      }
      if (length(p[p>0])==0) m <- "All probabilities are zero"
      if (m=="") {
        t <- cumsum(p/sum(p))
        t[t>0.99999999] <- 1.00000001
        x <- v[mapply(function(q,t) {which.max(cummax(q<t))},q,MoreArgs=list(t))]
      }
    } else if (s=="empi") {
      if (mode(v)!="numeric") {
        if (mode(try(scan(text=v,quiet=TRUE),silent=TRUE))=="numeric") {
          v <- c(scan(text=v,quiet=TRUE),recursive=TRUE) 
        } else {
          v <- c(scan(what=list(""),text=v,quiet=TRUE),recursive=TRUE) 
        }
      }  
      if (length(v)==0) v <- 0
      if (!is.na(lt)) lt <- min(v[which(v>=lt)])
      if (!is.na(ut)) ut <- max(v[which(v<=ut)])
      if (r=="y") {
        if (!is.na(lt)) v <- v[v>=lt] 
        if (!is.na(ut)) v <- v[v<=ut] 
      }
      if (length(v)==0) m <- "No empirical values"
      if (m=="") x <- v[round(0.5+length(v)*q)]
    } else if (s=="expo") {
      if (is.na(par1)) par1 <- 1
      if (is.na(par2)) par2 <- 0 
      if (!is.na(lt) && lt<par2) lt <- par2 
      if (!is.na(ut) && ut<par2) m <- "Upper expo. truncation below par2"
      if (!is.na(lt)&&!is.na(ut)&&lt>ut) m <- "Lower truncation above upper"
      if (par1<=0) { m <- paste("Invalid exponential parameter",par1)
      } else if (r=='y') { 
        if (!is.na(lt)) qlo <- pexp(lt-par2,par1)  else qlo <- 0
        if (!is.na(ut)) qhi <- pexp(ut-par2,par1)  else qhi <- 1 
        q <- qlo+q*(qhi-qlo)      
      }  
      if (m=="") x <- par2 + qexp(q,par1)
    } else if (s=="gamm") {
      if (is.na(par1)) par1 <- 1
      if (is.na(par2)) par2 <- 1
      if (is.na(par3)) par3 <- 0
      if (!is.na(lt) && lt<par3) lt <- par3
      if (!is.na(ut) && ut<par3) m <- "Upper gamma truncation below par3"
      if (!is.na(lt)&&!is.na(ut)&&lt>ut) m <- "Upper truncation above lower"
      if (par1<=0 | par2<=0) {
        m <- paste("Invalid gamma parameters",par1,par2,par3)
      } else if (r=='y') { 
        if (!is.na(lt)) qlo <- pgamma(lt-par3,par1,1/par2,par2) else qlo <- 0
        if (!is.na(ut)) qhi <- pgamma(ut-par3,par1,1/par2,par2) else qhi <- 1
        q <- qlo+q*(qhi-qlo)
      }  
      if (m=="") x <- par3 + qgamma(q,par1,1/par2,par2)
    } else if (s=="logn") {
      if (is.na(par1)) par1 <- 1
      if (is.na(par2)) par2 <- exp(1)
      if (is.na(par3)) par3 <- 0
      if (!is.na(lt) && lt<par3) lt <- par3
      if (!is.na(ut) && ut<par3) m <- "Upper lognormal truncation below par3"
      if (!is.na(lt)&&!is.na(ut)&&lt>ut) m <- "Upper truncation above lower" 
      if (par1<=0 | par2<1) {
        m <- paste("Invalid lognormal parameters",par1,par2,par3)
      } else if (r=='y' && par2>1) { 
        if (!is.na(lt)) qlo <- plnorm(lt-par3,log(par1),log(par2)) else qlo <- 0
        if (!is.na(ut)) qhi <- plnorm(ut-par3,log(par1),log(par2)) else qhi <- 1 
        q <- qlo+q*(qhi-qlo)      
      } 
      if (m=="") x <- par3 + qlnorm(q,log(par1),log(par2))
    } else if (s=="norm") {
      if (is.na(par1)) par1 <- 0
      if (is.na(par2)) par2 <- 1
      if (par2<0) {m <- paste("Invalid normal parameters",par1,par2)
      if (!is.na(lt)&&!is.na(ut)&&lt>ut) m <- "Lower truncation above upper"                 
      } else if (r=='y' && par2>0) { 
        if (!is.na(lt)) qlo <- pnorm(lt,par1,par2) else qlo <- 0
        if (!is.na(ut)) qhi <- pnorm(ut,par1,par2) else qhi <- 1 
        q <- qlo+q*(qhi-qlo)
      }  
      if (m=="") x <- qnorm(q,par1,par2)
    } else if (s=="poin") {
      if (is.na(par1)) par1 <- 0
      lt <- NA
      ut <- NA
      x <- rep(par1,length(q))
    } else if (s=="tria") {
      if (is.na(par1)&&is.na(par2)) { par1 <- 0; par2 <- 1
      } else if (is.na(par2)) { par2 <- par1+1
      } else if (is.na(par1))   par1 <- par2-1
      if (is.na(par3)) par3 <- (par1+par2)/2
      if (par3>par2) { t<-par2; par2<-par3; par3<-t }
      if (!is.na(lt) && lt>par2) m <- "Lower triangle truncation above par2"
      if (!is.na(ut) && ut<par1) m <- "Upper triangle truncation below par1"   
      if (!is.na(lt)&&!is.na(ut)&&lt>ut) m <- "Upper truncation above lower"
      if (par1>par2 | par3<par1 | par3>par2 ) {
        m <- paste("Invalid triangle parameters",par1,par2,par3)
      }  
      if (par1==par2) return(rep(par1,length(q)))
      p <- (par3-par1)/(par2-par1)
      if (r=='y') { 
        if (!is.na(lt) && lt>par1 && lt<=par2) { 
          if (lt==par3) { qlo <- p
          } else if (lt<par3) { qlo <-   (lt-par1)^2/((par2-par1)*(par3-par1))
          } else if (lt>par3) { qlo <- 1-(par2-lt)^2/((par2-par1)*(par2-par3))
          }                     
        } else qlo <- 0
        if (!is.na(ut) && ut>=par1 && ut<par2) {
          if (ut==par3) {qhi <- p
          } else if (ut<par3) { qhi <-   (ut-par1)^2/((par2-par1)*(par3-par1))
          } else if (ut>par3) { qhi <- 1-(par2-ut)^2/((par2-par1)*(par2-par3))
          }                     
        } else qhi <- 1
        q <- qlo+q*(qhi-qlo)
      }  
      if (m=="") { 
        x  <- par1 + sqrt(   q *(par2-par1)*(par3-par1))
        x2 <- par2 - sqrt((1-q)*(par2-par1)*(par2-par3))
        x[x2>par3] <- x2[x2>par3]
      }    
    } else if (s=="unif") {
      if (is.na(par1)) par1 <- 0
      if (is.na(par2)) par2 <- 1
      if (par2<par1) m <- paste("Invalid uniform parameters", par1,par2)
      if (!is.na(lt) && lt>par2) m <- "Lower uniform truncation above par2"
      if (!is.na(ut) && ut<par1) m <- "Upper uniform truncation below par1"
      if (!is.na(lt)&&!is.na(ut)&&lt>ut) m <- "Upper truncation above lower"    
      if (r=='y' && par2>par1) { 
        if (!is.na(lt)) qlo <- punif(lt,par1,par2) else qlo <- 0
        if (!is.na(ut)) qhi <- punif(ut,par1,par2) else qhi <- 1 
        q <- qlo+q*(qhi-qlo)
      }  
      if (m=="") x <- qunif(q,par1,par2)
    } else if (s=="weib") {
      if (is.na(par1)) par1 <- 1
      if (is.na(par2)) par2 <- 1
      if (is.na(par3)) par3 <- 0
      if (!is.na(lt) && lt<par3) lt <- par3
      if (!is.na(ut) && ut<par3) m <- "Upper Weibull truncation below par3"
      if (!is.na(lt)&&!is.na(ut)&&lt>ut) m <- "Upper truncation above lower"    
      if (par1<=0 | par2<=0) {
        m <- paste("Invalid Weibull parameters",par1,par2,par3) 
      }  
      if (r=='y') { 
        if (!is.na(lt)) qlo <- pweibull(lt-par3,par1,par2) else qlo <- 0
        if (!is.na(ut)) qhi <- pweibull(ut-par3,par1,par2) else qhi <- 1 
        q <- qlo+q*(qhi-qlo)
      } 
      if (m=="") x <- par3+qweibull(q,par1,par2)
    } else m <- paste("Unknown Distrib shape ",s)
    
    if (m != "") {cat("\n Error in distrib: ",m,"\n"); return(NULL)
    } else {
      if (!is.na(lt)) x <- mapply(max,lt,x)
      if (!is.na(ut)) x <- mapply(min,ut,x)   
      return(x)   
    } 
  }
  
  
  
  # eval.brand.list evaluates the products (brands) and formulations to be retained
  
  eval.brand.list = function(chem.list,puc.list,chem.fracs) {
    keep.puc <- rep(TRUE,length(puc.list))
    for (i in 1:length(puc.list)) {
      pchem <- unique(chem.fracs$dtxsid[chem.fracs$source.id==puc.list[i]])
      if (length(pchem %in% chem.list)==0) {
        cat ("\n PUC ",puc.list[i]," contains no modeled chemicals")
        keep.puc[i] <- FALSE
      }
    }  
    puc       <- puc.list[keep.puc]
    brands    <- vector("list",length(puc))
    for (i in 1:length(puc)) {
      brands[[i]] <- c(unique(chem.fracs$product_id[chem.fracs$source.id==puc[i]]))
    }
    brand.list  <- data.table(puc,brands)
    setorder(brand.list,puc)
    return(brand.list)
  }
  
  
  
  # eval.chem.list evaluates the chemicals to be retained
  
  eval.chem.list = function(chem.list,puc.list,chem.fracs) {
    keep.chem <- rep(TRUE,length(chem.list))
    chems     <- unique(chem.fracs$dtxsid)
    for (i in 1:length(chem.list)) {
      chem <- chem.list[i]
      if (!chem %in% chems) {
        cat("\n  Chem ",chem," not in any modeled PUC")
        keep.chem[i] <- FALSE
      } else {
        y <- chem.fracs[chem.fracs$dtxsid==chem]
        if (max(y$weight_fraction)==0) {
          cat ("\n  Chem ",chem," always zero")
          keep.chem[i] <- FALSE
        }  
      }      
    }  
    chem.list <- unique(chem.list[keep.chem])
    return(chem.list)
  }
  
  
  
  # eval.chem.props returns chemical properties sampled from distributions
  
  eval.chem.props = function(fug.cvars,chem.list,ran.vars,q) {
    nc         <- length(chem.list)
    molwt      <- vector("numeric",nc)
    vapor      <- vector("numeric",nc)
    solub      <- vector("numeric",nc)
    kow        <- vector("numeric",nc)
    decay.air  <- vector("numeric",nc)
    decay.sur  <- vector("numeric",nc)
    diffus.air <- vector("numeric",nc)
    conc.out   <- vector("numeric",nc)
    qf         <- as.data.frame(q)[ran.vars[[3]]]
    qx         <- data.table(matrix(qf,nrow=nc,ncol=ncol(qf)/nc,byrow=TRUE))
    setnames(qx,unique(unlist(lapply(names(qf),strip.n))))
    if(g$save.r.objects=="y") write.csv(unlist(qx),paste0(g$out,"/Temp/qx_",run.name,".csv"))
    for (c in 1:nc) {
      cvars <- fug.cvars[fug.cvars$dtxsid==chem.list[c]]
      if (nrow(cvars)>0) {
        qc            <- lapply(qx[c],as.numeric)
        molwt[c]      <- cvars$molwt
        vaporGM       <- cvars$vapor
        if(is.null(vaporGM)) vaporGM <- 1
        vapor[c]      <- distrib("logn",par1=vaporGM,par2=2,lt=vaporGM/100,ut=min(vaporGM*100,1E5),q=qc$vapor)
        solubGM       <- cvars$solub
        if(is.null(solubGM))    solubGM <- 1
        solub[c]      <- distrib("logn",par1=solubGM,par2=1.3,lt=solubGM/5,ut=solubGM*5,q=qc$solub)
        kowGM         <- cvars$kow
        if(is.null(kowGM))      kowGM <- 10
        kow[c]        <- distrib("logn",par1=kowGM,par2=1.5,lt=kowGM/10,ut=kowGM*10,q=qc$kow)
        dc.airGM      <- cvars$decay.a
        if(is.null(dc.airGM))   dc.airGM <- 0.01
        decay.air[c]  <- distrib("logn",par1=dc.airGM,par2=1.5,lt=dc.airGM/10,ut=dc.airGM*10,q=qc$decay.air)
        dc.surGM      <- cvars$decay.s
        if(is.null(dc.surGM))   dc.surGM <- 0.001
        decay.sur[c]  <- distrib("logn",par1=dc.surGM,par2=1.5,lt=dc.surGM/10,ut=dc.surGM*10,q=qc$decay.sur)
        diffusGM      <- (2.05*(1/29+1/molwt[c])^0.5)/molwt[c]^0.33 * 86400 / 10000 
        diffus.air[c] <- distrib("logn",par1=diffusGM,par2=1.3,lt=diffusGM/5,ut=diffusGM*5,q=qc$diffus.air)
        backGM        <- cvars$backgm
        backGSD       <- cvars$backgsd
        if (backGM>0) {
          conc.out[c] <- distrib("logn",par1=backGM,par2=backGSD,lt=backGM/10,ut=backGM*10,q=qc$conc.out)
        } else conc.out[c] <- 0
      } else {
        molwt[c]      <- 100
        vapor[c]      <- 1
        solub[c]      <- 1
        kow[c]        <- 10
        decay.air[c]  <- 0.01
        decay.sur[c]  <- 0.001
        diffus.air[c] <- 0.01
        conc.out[c]   <- 0
      }  
    }
    props <- data.table(chem.list,molwt,vapor,solub,kow,decay.air,decay.sur,diffus.air,conc.out)
    setnames(props,names(props),c("dtxsid","molwt","vapor","solub","kow","decay.air","decay.sur","diffus.air","conc.out"))
    if(g$prog=="y") cat("\n  Evaluating chem props complete")
    return(props)
  }
  
  
  
  # eval.chem.release produces an hourly calendar with chemical releases inside house
  
  eval.chem.release = function(pucs,use.data,use.chem,compart.list,diary,nc,calendar.hours) {
    np        <- nrow(pucs)
    air       <- matrix(use.data[,3,compart.list=="fia",],nrow=np,ncol=nc)
    sur       <- matrix(use.data[,3,compart.list=="fis",],nrow=np,ncol=nc)
    tot       <- as.data.table(use.chem)
    setnames(tot,str_c(rep("tot",nc),rep(1:nc)))
    puc.data  <- data.table(pucs$source.id,air,sur,tot)
    setnames(puc.data,1:(1+2*nc),c("source.id",str_c(rep("air",nc),rep(1:nc)),str_c(rep("sur",nc),rep(1:nc))))
    cair      <- as.data.table(matrix(0,nrow=8736,ncol=nc))
    csur      <- as.data.table(matrix(0,nrow=8736,ncol=nc))
    ctot      <- as.data.table(matrix(0,nrow=8736,ncol=nc))
    rel.ind   <- rep(0,8736)
    rel.tot   <- rep(0,8736)
    setnames(cair,str_c(rep("air",nc),rep(1:nc)))
    setnames(csur,str_c(rep("sur",nc),rep(1:nc)))
    setnames(ctot,str_c(rep("tot",nc),rep(1:nc)))
    d         <- select(diary[diary$mass!=0],source.id,row,person,daynum,hour,mass)
    d[!d$source.id %in% pucs$source.id]$source.id <- "none"
    d$hournum <- as.integer(24*(d$daynum-1)+d$hour)
    setkey(d,source.id)#sorting by source.id, in ascending order
    setkey(puc.data,source.id)
    dpuc      <- as.data.frame(inner_join(puc.data,d,by="source.id"))
    hournum   <- dpuc$hournum
    yair      <- select_vars(names(puc.data),starts_with("air"))
    ysur      <- select_vars(names(puc.data),starts_with("sur"))
    ytot      <- select_vars(names(puc.data),starts_with("tot"))
    pair      <- as.data.table(dpuc[yair])
    psur      <- as.data.table(dpuc[ysur])
    ptot      <- as.data.table(dpuc[ytot])
    hours     <- unique(hournum)
    for (i in 1:length(hours)) {
      hr <- hours[i]
      set(cair,i=hr,1:nc,lapply(pair[hournum==hr],sum))
      set(csur,i=hr,1:nc,lapply(psur[hournum==hr],sum))
      set(ctot,i=hr,1:nc,lapply(ptot[hournum==hr],sum))
      rel.ind[hr] <- sum(cair[hr],csur[hr])
      rel.tot[hr] <- sum(ctot[hr])
    } 
    chem.release <- as.data.table(data.frame(calendar.hours,rel.ind,rel.tot,cair,csur,ctot))
    if (g$prog=="y") cat("\n  Evaluating chemical release complete")
    return(chem.release)   
  }
  
  
  
  # eval.chem.totals sets the annual household usage of each chemical in mg
  
  eval.chem.totals = function(chem.release,chem.list,house.num) {
    if (!is.null(chem.release)) {
      x <- select(as.data.frame(chem.release),starts_with("tot"))
      y <- as.data.table(matrix(colSums(x),nrow=1,ncol=length(chem.list)))
    } else y <- as.data.table(matrix(0,nrow=1,ncol=length(chem.list)))
    setnames(y,chem.list)
    write.csv(y,paste0(g$out,"/Chem/House_",house.num,".csv"),row.names=TRUE)
    return(y)
  }
  
  
  
  
  # eval.dermal.rates determines the rate constants for removal processes from skin 
  
  eval.dermal.rates = function(fug.cvars,hp,prime,pucs.areas,prod.chem,nc,house.num) {
    np    <- nrow(pucs.areas)                                             # last one is for indirect
    kaw   <- 1/(fug.cvars$solub/fug.cvars$vapor*8.314*hp$temp)            # inverse of SHEDS convention
    phi.a <- 8 * kaw
    phi.w <- 20.6/fug.cvars$molwt^0.4757
    kp    <- fug.cvars$kp /100                                            # convert cm/hr to m/hr    
    kh    <- prime$hand.mouth
    kr    <- 0         # other removal is zero for now
    # hw is number per day, all others are rates per hour
    hw    <- prime$hand.wash
    # ka and ks require product layer thickness on skin
    # first calulate affected skin area by PUC (in cm2 from pophouse file)
    #if (prime$age<18)  prod.area <- pucs.areas$f.child * prime$skin.area
    #if (prime$age>=18) prod.area <- pucs.areas$f.adult * prime$skin.area
    #prod.area[is.na(prod.area)]  <- 1
    #prod.area                    <- pmax(prod.area, 1)                  # set minimum of 1 cm2  
    if (prime$age<18)  frac.area <- pmax(pucs.areas$f.child,0.01)        # set minimum to 1% of total
    if (prime$age>=18) frac.area <- pmax(pucs.areas$f.adult,0.01)
    frac.area[is.na(frac.area)]  <- 0.01
    prod.area                    <- frac.area * prime$skin.area
    # volume in cm3 = product mass in grams * fraction applied to skin (assume density=1 g/cm3)
    prod.vol  <- pucs.areas$mass * c(prod.chem$fsk,1)
    # thickness is volume / affected skin area (cm3/cm2)
    prod.thick <- prod.vol / prod.area 
    h            <- prod.thick / 100        # convert cm to meters
    h[length(h)] <- pucs.areas$h[length(h)]
    h            <- pmax(1E-5,h)            # minimum of 10 microns 
    source.id <- rep(pucs.areas$source.id,nc)
    dtxsid    <- unlist(lapply(fug.cvars$dtxsid,rep,np))
    rates     <- data.table(source.id,dtxsid,matrix(0,length(dtxsid),ncol=6))
    setnames(rates,c("source.id","dtxsid","hw","kh","kr","ka","ks","h"))
    rates$hw <- hw
    rates$kh <- kh
    rates$kr <- kr
    for (c in 1:length(kaw)) {
      for (p in 1:np) {
        rates$ka[p+(c-1)*np] <- 1/(h[[p]]/phi.a[[c]] + h[[p]]/phi.w[[c]])  # rate into air
        rates$ks[p+(c-1)*np] <- 1/(h[[p]]/kp[[c]]    + h[[p]]/phi.w[[c]])  # rate into skin
      }
    }
    rates$h  <- h
    rates[is.na(rates)] <- 0
    if (g$prog=="y") cat("\n  Evaluating dermal rates complete")
    if (g$save.r.objects=="y") write.csv(rates,paste0(g$out,"/Temp/dermal_rates_",house.num,"_",run.name,".csv"))
    return(rates)
  }
  
  
  
  
  # eval.direct produces a daily summary of direct exposure variables
  
  eval.direct = function(d,use.data,use.chem,pucs,dermal.rates,compart.list,prime,puc.wipe.rinse,chem.list,fug.cvars,chem.totals)  {
    chems        <- as.data.frame(chem.totals)
    # First, evaluate the exposure variables for each combination of product and chemical
    skin1        <- use.data[,1,compart.list=="fsk",]
    skin2        <- use.data[,2,compart.list=="fsk",]
    skin3        <- use.data[,3,compart.list=="fsk",]
    inair1       <- use.data[,1,compart.list=="fia",]
    inair2       <- use.data[,2,compart.list=="fia",]
    gut          <- use.data[,3,compart.list=="fgi",]
    nc           <- length(chem.list)
    np           <- nrow(pucs)
    use.avg.derm <- as.data.frame(matrix((skin1+skin2)/2,nrow=np,ncol=nc))
    use.avg.hand <- use.avg.derm * puc.wipe.rinse$fhands
    use.avg.body <- use.avg.derm * puc.wipe.rinse$fbody
    use.avg.air  <- as.data.frame(matrix((inair1+inair2)/2,nrow=np,ncol=nc))
    use.gut      <- as.data.frame(matrix(gut,nrow=np,ncol=nc))
    post.derm    <- as.data.frame(matrix(skin3,nrow=np,ncol=nc))
    post.hand    <- post.derm * puc.wipe.rinse$fhands
    post.body    <- post.derm * puc.wipe.rinse$fbody
    vol.cloud    <- 2                                        # personal cloud size
    aer.cloud    <- 10                                       # cloud exchanges air 10 times per hour
    dur.cloud    <- pucs$hand.dur/60                         # handling time in hours
    rates        <- as.data.table(dermal.rates)[!dermal.rates$source.id=="Indirect"]
    f.hands      <- puc.wipe.rinse$fhands
    f.hands[is.na(f.hands)] <- 0.05
    hand.area    <- 0.05*prime$skin.area         
    body.area    <- 0.95*prime$skin.area
    vars <- c("f.use.lost","use.avg.air","use.avg.derm","use.hands.exp","use.body.exp","use.hands.air",
              "use.body.air","use.body.abs","use.hands.load","use.body.load","use.dermal.abs","use.air.mass",
              "use.air.conc","use.inhal.exp","use.inhal.mass","use.hands.abs",
              "use.inhal.abs","use.inges.exp","use.inges.abs",
              "f.hands.lost","f.body.lost","post.dermal.exp","post.hands.exp","post.body.exp","post.hands.abs",
              "post.hands.air","post.inges.hand","post.body.air","post.body.abs",
              "post.hands.load","post.body.load","post.dermal.abs",
              "post.air.rate","post.air.conc","post.inhal.exp","post.inhal.mass",
              "post.inhal.abs","post.inges.abs")
    zero   <- data.table(matrix(0,nrow=nrow(pucs),ncol=length(vars)))
    setnames(zero,vars)
    dp     <- data.table(pucs$source.id,f.hands,pucs$hand.dur,pucs$hand.dur/60,vol.cloud,aer.cloud)
    setnames(dp,c("source.id","f.hands","use.dur.min","use.dur.hr","vol.cloud","aer.cloud"))
    if (g$prog=="y") cat("\n     direct initiation complete...")
    for (c in 1:nc) {
      dp <- cbind(dp,zero)
      if (chems[c]>0) {
        if (g$prog=="y") cat(paste0("\n   Chemical ",c,"  ",chem.list[c]))
        ka <- pmax(0.00001,rates$ka[rates$dtxsid==chem.list[c]])
        ks <- pmax(0.00001,rates$ks[rates$dtxsid==chem.list[c]])
        kh <- pmax(0.00001,rates$kh[rates$dtxsid==chem.list[c]])
        kr <- pmax(0.00001,rates$kr[rates$dtxsid==chem.list[c]])
        hw <- pmax(0.00001,rates$hw[rates$dtxsid==chem.list[c]])
        dp$f.use.lost       <- 1-exp(-(ka+ks)*dp$use.dur.hr)    # fraction of hand loading lost during product usage time
        dp$use.avg.air      <- use.avg.air[c]
        dp$use.avg.derm     <- use.avg.derm[c]
        dp$use.hands.exp    <- use.avg.derm[c] * dp$f.hands
        dp$use.body.exp     <- use.avg.derm[c] * (1-dp$f.hands)
        dp$use.hands.air    <- dp$use.hands.exp * ka/(ka+ks)* dp$f.use.lost
        dp$use.hands.abs    <- dp$use.hands.exp * ks/(ka+ks)* dp$f.use.lost
        dp$use.body.air     <- dp$use.body.exp  * ka/(ka+ks)* dp$f.use.lost
        dp$use.body.abs     <- dp$use.body.exp  * ks/(ka+ks)* dp$f.use.lost
        dp$use.hands.load   <- dp$use.hands.exp / hand.area
        dp$use.body.load    <- dp$use.body.exp  / body.area
        dp$use.dermal.abs   <- dp$use.hands.abs + dp$use.body.abs
        dp$use.air.mass     <- dp$use.avg.air + dp$use.hands.air + dp$use.body.air
        dp$use.air.conc     <- dp$use.air.mass/(vol.cloud*aer.cloud*dur.cloud)
        dp$use.inhal.exp    <- dp$use.air.conc * dp$use.dur.hr/24             # converted to daily avg
        dp$use.inhal.mass   <- dp$use.inhal.exp * prime$basal.vent * pucs$met
        dp$use.inges.exp    <- use.gut[c]                                     # amount in fgi compartment
        # for spray products, reassign 75% of inhaled mass to ingestion
        dp$use.inges.exp[pucs$spray]  <- 0.75*dp$use.inhal.mass[pucs$spray] + dp$use.inges.exp[pucs$spray]  
        dp$use.inhal.mass[pucs$spray] <- 0.25*dp$use.inhal.mass[pucs$spray] 
        dp$use.inges.abs    <- dp$use.inges.exp * fug.cvars$fabs[c]
        dp$use.inhal.abs    <- dp$use.inhal.mass * 0.16                       # assume 16% absorption 
        
        dp$f.hands.lost     <- 1-exp(-(ka+ks+kh+kr)*8/hw)                     # dur on hands is (8/handwashes) hrs
        dp$f.body.lost      <- 1-exp(-(ka+ks+kr)*8)                           # dur on body is 8 hrs
        dp$post.dermal.exp  <- post.derm[c]
        dp$post.hands.exp   <- post.derm[c] * dp$f.hands
        dp$post.body.exp    <- post.derm[c] * (1-dp$f.hands)
        dp$post.hands.abs   <- dp$post.hands.exp * ks/(ka+ks+kh+kr) * dp$f.hands.lost
        dp$post.hands.air   <- dp$post.hands.exp * ka/(ka+ks+kh+kr) * dp$f.hands.lost
        dp$post.inges.hand  <- dp$post.hands.exp * kh/(ka+ks+kh+kr) * dp$f.hands.lost
        dp$post.body.air    <- dp$post.body.exp  * ka/(ka+ks+kr)    * dp$f.body.lost
        dp$post.body.abs    <- dp$post.body.exp  * ks/(ka+ks+kr)    * dp$f.body.lost
        dp$post.hands.load  <- dp$post.hands.exp / hand.area
        dp$post.body.load   <- dp$post.body.exp  / body.area
        dp$post.dermal.abs  <- dp$post.hands.abs+dp$post.body.abs
        dp$post.air.rate    <- (dp$post.hands.air + dp$post.body.air)/8       # hourly emission rate into cloud
        dp$post.air.conc    <- dp$post.air.rate/(vol.cloud*aer.cloud)         # inflow in 1 hour / ouflow in 1 hour 
        dp$post.inhal.exp   <- dp$post.air.conc * 8/24                        # convert to daily avg
        dp$post.inhal.mass  <- dp$post.inhal.exp * prime$basal.vent * 2.2     # assume post-use mean MET is 2.2
        dp$post.inhal.abs   <- dp$post.inhal.mass * 0.16                      # assume 16% absorption
        dp$post.inges.abs   <- dp$post.inges.hand * fug.cvars$fabs[c]
      }
      setnames(dp,vars,str_c(vars,c)) 
    }
    dchem <- as.data.table(data.frame(dp,use.chem))
    setkey(dchem,source.id)
    
    # Now loop over the diary events, applying the above results as products are used
    d <- d[d$source.id %in% pucs$source.id]
    setkey(d,source.id)
    ddata     <- left_join(d,dchem,by="source.id")
    derm.exp  <- matrix(0,nrow=364,ncol=nc)
    derm.max  <- matrix(0,nrow=364,ncol=nc)
    derm.abs  <- matrix(0,nrow=364,ncol=nc)
    inhal.exp <- matrix(0,nrow=364,ncol=nc)
    inhal.mass<- matrix(0,nrow=364,ncol=nc)
    inhal.max <- matrix(0,nrow=364,ncol=nc)
    inhal.abs <- matrix(0,nrow=364,ncol=nc)
    inges.exp <- matrix(0,nrow=364,ncol=nc)
    inges.abs <- matrix(0,nrow=364,ncol=nc)
    release   <- matrix(0,nrow=364,ncol=nc)
    nam.de    <- str_c(rep("dir.derm.exp",nc),1:nc)
    nam.dm    <- str_c(rep("dir.derm.max",nc),1:nc)
    nam.da    <- str_c(rep("dir.derm.abs",nc),1:nc)
    nam.ihe   <- str_c(rep("dir.inhal.exp",nc),1:nc)
    nam.ihmas <- str_c(rep("dir.inhal.mass",nc),1:nc)
    nam.ihm   <- str_c(rep("dir.inhal.max",nc),1:nc)
    nam.iha   <- str_c(rep("dir.inhal.abs",nc),1:nc)
    nam.ige   <- str_c(rep("dir.ingest.exp",nc),1:nc)
    nam.iga   <- str_c(rep("dir.ingest.abs",nc),1:nc)
    nam.rel   <- str_c(rep("dir.release",nc),1:nc)
    uhe <- select_vars(names(ddata),starts_with("use.hands.exp"))
    ube <- select_vars(names(ddata),starts_with("use.body.exp"))
    phe <- select_vars(names(ddata),starts_with("post.hands.exp"))
    pbe <- select_vars(names(ddata),starts_with("post.body.exp"))
    uha <- select_vars(names(ddata),starts_with("use.hands.abs"))
    uba <- select_vars(names(ddata),starts_with("use.body.abs"))
    pha <- select_vars(names(ddata),starts_with("post.hands.abs"))
    pba <- select_vars(names(ddata),starts_with("post.body.abs"))
    uhl <- select_vars(names(ddata),starts_with("use.hands.load"))
    ubl <- select_vars(names(ddata),starts_with("use.body.load"))
    phl <- select_vars(names(ddata),starts_with("post.hands.load"))
    pbl <- select_vars(names(ddata),starts_with("post.body.load"))
    uie <- select_vars(names(ddata),starts_with("use.inhal.exp"))
    pie <- select_vars(names(ddata),starts_with("post.inhal.exp"))
    uim <-  select_vars(names(ddata),starts_with("use.inhal.mass"))
    pim <- select_vars(names(ddata),starts_with("post.inhal.mass"))
    uac <- select_vars(names(ddata),starts_with("use.air.conc"))
    pac <- select_vars(names(ddata),starts_with("post.air.conc"))
    uia <- select_vars(names(ddata),starts_with("use.inhal.abs"))
    pia <- select_vars(names(ddata),starts_with("post.inhal.abs"))
    pgh <- select_vars(names(ddata),starts_with("post.inges.hand"))
    ugm <- select_vars(names(ddata),starts_with("use.inges.mass"))
    uga <- select_vars(names(ddata),starts_with("use.inges.abs"))
    pga <- select_vars(names(ddata),starts_with("post.inges.abs"))
    tot <- select_vars(names(ddata),starts_with("tot"))
    if (g$prog=="y") cat("\n     direct starting loop over chemicals and days..")
    ndaysused <- nrow(ddata)
    if (ndaysused>0) {
      for (i in 1:ndaysused) {
        dd  <- as.data.frame(ddata[i])
        day <- d$daynum[i]
        for (c in 1:nc) {
          if (chems[c]>0) {
            derm.exp[day,c]   <- sum(dd[uhe][c],dd[ube][c],dd[phe][c],dd[pbe][c])
            derm.max[day,c]   <- max(0,unlist(dd[uhl][c]),unlist(dd[ubl][c]),unlist(dd[phl][c]),unlist(dd[pbl][c]))
            derm.abs[day,c]   <- sum(dd[uha][c],dd[uba][c],dd[pha][c],dd[pba][c])
            inhal.exp[day,c]  <- sum(dd[uie][c],dd[pie][c])
            inhal.mass[day,c] <- sum(dd[uim][c],dd[pim][c])
            inhal.max[day,c]  <- max(0,unlist(dd[uac][c]),unlist(dd[pac][c]),unlist(dd[uie][c]),unlist(dd[pie][c]))
            inhal.abs[day,c]  <- sum(dd[uia][c],dd[pia][c])
            inges.exp[day,c]  <- sum(dd[pgh][c])
            inges.abs[day,c]  <- sum(dd[uga][c],dd[pga][c])
            release[day,c]    <- sum(dd[tot][c])
          }
        }  
      }
    }
    direct <- as.data.table(data.frame(1:364,derm.exp,derm.max,derm.abs,inhal.exp,inhal.mass,inhal.max,inhal.abs,inges.exp,inges.abs,release))
    setnames(direct,c("daynum",nam.de,nam.dm,nam.da,nam.ihe,nam.ihmas,nam.ihm,nam.iha,nam.ige,nam.iga,nam.rel))
    if (g$prog=="y") cat("\n  Evaluating direct exposures complete")
    return(direct)
  }
  
  
  
  # eval.env.impact produces a daily summary of emissions from the house
  
  eval.env.impact = function(diary,use.data,pucs,fug.day,compart.list,chem.totals,nc) {
    if (!is.null(use.data)) {
      chems   <- as.data.frame(chem.totals)
      d       <- diary[diary$source.id %in% pucs$source.id]
      np      <- nrow(pucs)
      out.sur <- matrix(use.data[,3,compart.list=="fos",],nrow=np,ncol=nc)
      out.air <- matrix(use.data[,3,compart.list=="foa",],nrow=np,ncol=nc)
      drain   <- matrix(use.data[,3,compart.list=="fdr",],nrow=np,ncol=nc)
      waste   <- matrix(use.data[,3,compart.list=="fws",],nrow=np,ncol=nc)
      rels    <- data.table(pucs$source.id,drain,waste,out.sur,out.air)
      nam.os  <- str_c(rep("out.sur",nc),1:nc)
      nam.oa  <- str_c(rep("out.air",nc),1:nc)
      nam.dr  <- str_c(rep("drain",nc),1:nc)
      nam.ws  <- str_c(rep("waste",nc),1:nc)
      setnames(rels,c("source.id",nam.dr,nam.ws,nam.os,nam.oa))
      setkey(d,source.id)
      setkey(rels,source.id)
      drel <- left_join(d,rels,by="source.id")
      dir.os <- matrix(0,nrow=364,ncol=nc)
      dir.oa <- matrix(0,nrow=364,ncol=nc)
      dir.dr <- matrix(0,nrow=364,ncol=nc)
      dir.ws <- matrix(0,nrow=364,ncol=nc)
      os  <- select_vars(names(drel),starts_with("out.sur"))
      oa  <- select_vars(names(drel),starts_with("out.air"))
      dr  <- select_vars(names(drel),starts_with("drain"))
      ws  <- select_vars(names(drel),starts_with("waste"))
      fwn <- select_vars(names(fug.day),starts_with("win"))
      fws <- select_vars(names(fug.day),starts_with("was"))
      for (c in 1:nc) {
        if(chems[c]>0) {
          for (i in 1:364) {
            dd          <- as.data.frame(drel[drel$daynum==i])
            if (nrow(dd)>0) {
              fd          <- as.data.frame(fug.day[fug.day$daynum==i])
              dir.os[i,c] <- sum(dd[os][c])
              dir.oa[i,c] <- sum(dd[oa][c],fd[fwn][c])
              dir.dr[i,c] <- sum(dd[dr][c])
              dir.ws[i,c] <- sum(dd[ws][c],fd[fws][c])
            }  
          }  
        }
      }
      env.impact <- data.table(1:364,dir.os,dir.oa,dir.dr,dir.ws)
      setnames(env.impact,c("daynum",os,oa,dr,ws))
    } else {
      env.impact <- data.table(0)
    }  
    if (g$prog=="y") cat("\n  Evaluating environmental impacts complete")
    return(env.impact)
  }
  
  
  
  
  # eval.flows detetmines the flow rate constants between house compartments
  
  eval.flows = function(hp,cp) { 
    cp$vol.air      <- hp$area.sur * hp$height                            # vol.air [m3]
    cp$aer.out      <- hp$aer.out                                         # aer.out [1/day]
    cp$sm.mass.air  <- cp$vol.air * hp$sm.load.air                        # sm.mass.air [ug]
    cp$lg.mass.air  <- cp$vol.air * hp$lg.load.air                        # lg.mass.air [ug]
    cp$sm.mass.sur  <- hp$area.sur * hp$sm.load.sur                       # sm.mass.sur [ug]
    cp$lg.mass.sur  <- hp$area.sur * hp$lg.load.sur                       # lg.mass.sur [ug]
    cp$sm.clean.sur <- max(hp$sm.clean.sur,hp$sm.depos*hp$sm.load.air/hp$sm.load.sur-hp$sm.resus)
    cp$lg.clean.sur <- max(hp$lg.clean.sur,hp$lg.depos*hp$lg.load.air/hp$lg.load.sur-hp$lg.resus)
    # clean.sur [1/day], depos [m/day], load.air [ug/m3], load.sur [ug/m2], resus [1/day]
    cp$sm.clean.air <- pmax(hp$sm.clean.air,(hp$sm.resus*cp$sm.mass.sur-hp$sm.depos*hp$sm.load.air*hp$area.sur)/cp$sm.mass.air)
    cp$lg.clean.air <- pmax(hp$lg.clean.air,(hp$lg.resus*cp$lg.mass.sur-hp$lg.depos*hp$lg.load.air*hp$area.sur)/cp$lg.mass.air)
    # clean.air [1/day], resus [1/day], mass.sur [ug], depos [m/day], load.air [ug/m3], area.sur [m2], mass.air [ug]
    cp$ug.mol       <- 1E6*cp$molwt                                       # ug.mol [ug/mol]
    cp$z.air        <- 1/(8.314*hp$temp)                                  # z.air [mol/(Pa m3)]
    cp$zvb.air      <- cp$z.air * cp$vol.air * cp$ug.mol                  # zvb.air [ug/Pa]
    cp$z.sur        <- cp$z.air * 82500 / (cp$vapor^0.65)                 # z.sur [mol/(Pa m3)]
    cp$zvb.sur      <- cp$z.sur * hp$area.sur * hp$thick.sur * cp$ug.mol  # zvb.sur [ug/Pa]
    cp$sm.kp        <- 1.662E-12 * cp$kow * hp$sm.carb.f * cp$solub / (cp$vapor * cp$z.air)
    cp$lg.kp        <- 1.662E-12 * cp$kow * hp$lg.carb.f * cp$solub / (cp$vapor * cp$z.air)
    # kp [m3/ug], 1.662E-12 [m3/ug], kow [-], carb.f [-], solub [mol/m3], vapor [Pa], z.air [mol/(Pa m3)]
    cp$sm.zv.air    <- cp$zvb.air * cp$sm.kp * hp$sm.load.air             # sm.zv.air [ug/Pa]
    cp$lg.zv.air    <- cp$zvb.air * cp$lg.kp * hp$lg.load.air             # lg.zv.air [ug/Pa]
    cp$sm.cap       <- cp$z.air * cp$sm.kp * cp$ug.mol                    # sm.cap [1/Pa]
    cp$lg.cap       <- cp$z.air * cp$lg.kp * cp$ug.mol                    # lg.cap [1/Pa]
    cp$sm.zv.sur    <- cp$sm.cap * cp$sm.mass.sur                         # sm.zv.sur [ug/Pa]
    cp$lg.zv.sur    <- cp$lg.cap * cp$lg.mass.sur                         # lg.zv.sur [ug/Pa]
    cp$zv.air       <- cp$zvb.air + cp$sm.zv.air + cp$lg.zv.air           # zv.air [ug/Pa]
    cp$zv.sur       <- cp$zvb.sur + cp$sm.zv.sur + cp$lg.zv.sur           # zv.sur [ug/Pa]
    cp$izv.air      <- pmin(1E100,1/cp$zv.air)                            # izv.air [Pa/ug]
    cp$izv.sur      <- pmin(1E100,1/cp$zv.sur)                            # izv.sur [Pa/ug]
    cp$yaf          <- pmin(cp$diffus.air*cp$z.air/hp$thick.bou,0.0135/(cp$vapor^0.32))        # yaf [mol/(m2-Pa-day)]
    cp$cln.air      <- cp$izv.air*(cp$sm.zv.air*cp$sm.clean.air+cp$lg.zv.air*cp$lg.clean.air)  # cln.air [1/day]
    cp$cln.sur      <- cp$izv.sur*(cp$sm.zv.sur*cp$sm.clean.sur+cp$lg.zv.sur*cp$lg.clean.sur)  # cln.sur [1/day]
    cp$sm.dep       <- cp$izv.air * hp$area.sur * hp$sm.load.air * hp$sm.depos * cp$sm.cap     # sm.dep [1/day]
    cp$lg.dep       <- cp$izv.air * hp$area.sur * hp$lg.load.air * hp$lg.depos * cp$lg.cap     # lg.dep [1/day]
    cp$dep          <- cp$sm.dep + cp$lg.dep                                                   # dep [1/day]
    cp$res          <- cp$izv.sur*(cp$sm.mass.sur*hp$sm.resus+cp$lg.mass.sur*hp$lg.resus)      # res [1/day]
    cp$diff.air     <- cp$izv.air * cp$ug.mol * hp$area.sur * cp$yaf                           # diff.air [1/day]
    cp$diff.sur     <- cp$izv.sur * cp$ug.mol * hp$area.sur * cp$yaf                           # diff.sur [1/day]
    cp$a <- hp$aer.out + cp$decay.air + cp$cln.air + cp$dep + cp$diff.air                      # a [1/day]
    cp$b <- cp$res + cp$diff.sur                                                               # b [1/day]
    cp$c <- cp$dep + cp$diff.air                                                               # c [1/day]
    cp$d <- cp$decay.sur + cp$cln.sur + cp$res + cp$diff.sur                                   # d [1/day]
    if (g$prog=="y") cat("\n  Evaluating flows complete")
    return(cp)
  }
  
  
  
  
  # eval.fug.concs.an solves the fugacity equations analytically
  
  eval.fug.concs.an = function(chem.release,flows,indoor.hrs,indoor.gaps) {
    if (g$prog=="y") cat("\n  Starting fugacity calculations")
    nc          <- nrow(flows)
    air.mass    <- as.data.table(matrix(0,nrow=8736,ncol=nc))
    sur.mass    <- as.data.table(matrix(0,nrow=8736,ncol=nc))
    waste.mass  <- as.data.table(matrix(0,nrow=364,ncol=nc))
    window.mass <- as.data.table(matrix(0,nrow=364,ncol=nc))
    air <- rep(0,nc) 
    sur <- rep(0,nc)
    was <- rep(0,nc)
    win <- rep(0,nc)
    if (!is.null(chem.release)) {
      yair        <- select_vars(names(chem.release),starts_with("air"))
      ysur        <- select_vars(names(chem.release),starts_with("sur"))
      df          <- as.data.frame(chem.release)
      new.air     <- as.data.table(df[yair])                # chemical mass in mg 
      new.sur     <- as.data.table(df[ysur])                # chemical mass in mg
    } else {
      new.air <- rep(0,8736)
      new.sur <- rep(0,8736)
    }
    air <- rep(0,nc) 
    sur <- rep(0,nc)
    was <- rep(0,nc)
    win <- rep(0,nc)
    fac   <- 24
    fromair <- flows$a/fac                               # get hourly flow rates
    toair   <- flows$b/fac
    tosur   <- flows$c/fac
    fromsur <- flows$d/fac
    rterm   <- sqrt(fromair^2+4*toair*tosur-2*fromair*fromsur+fromsur^2)
    lam1    <- (fromair+fromsur+rterm)/2                 # eigenvalues of J
    lam2    <- (fromair+fromsur-rterm)/2                 # eigenvalues of J
    maxlag  <- max(8735,indoor.gaps)+1
    hrvec   <- 1:maxlag
    elt1    <- as.data.frame(matrix(0,nrow=maxlag,ncol=nc))
    elt2    <- as.data.frame(matrix(0,nrow=maxlag,ncol=nc))
    for (c in 1:nc) {
      elt1[c] <- as.data.frame.vector(lapply(-lam1[c]*hrvec,exp))
      elt2[c] <- as.data.frame.vector(lapply(-lam2[c]*hrvec,exp))
    }
    towin   <- flows$aer.out/fac
    cl.air  <- flows$cln.air/fac
    cl.sur  <- flows$cln.sur/fac
    # new code to handle outdoor background as a source
    cout    <- flows$conc.out/1000                                     # convert from ug/m3 to mg/m3     
    src.air <- cout*flows$aer.out*flows$vol.air                        # units are mg/day
    mc.air  <- src.air*flows$d / (flows$a*flows$d - flows$b*flows$c)   # units are mg
    mc.sur  <- src.air*flows$c / (flows$a*flows$d - flows$b*flows$c)   # units are mg
    if(length(indoor.hrs)>0) {
      for (hrnum in 1:length(indoor.hrs)) {
        hr  <- as.integer(indoor.hrs[hrnum])
        hr1 <- as.integer(max(1,hr-1))
        ma0 <- air.mass[hr1] + new.air[hr]
        ms0 <- sur.mass[hr1] + new.sur[hr]
        k1a <- (2*tosur*ma0+(fromair-fromsur-rterm)*ms0)*(rterm+fromair-fromsur)/(4*tosur*rterm)
        k2a <- (2*tosur*ma0+(fromair-fromsur+rterm)*ms0)*(rterm-fromair+fromsur)/(4*tosur*rterm)
        k1s <- -(2*tosur*ma0+(fromair-fromsur-rterm)*ms0)/(2*rterm)
        k2s <-  (2*tosur*ma0+(fromair-fromsur+rterm)*ms0)/(2*rterm)
        n  <- indoor.gaps[hrnum]+1
        if (n>0) {
          for (lag in 1:n){
            t1 <- as.numeric(unlist(elt1[lag,]))
            t2 <- as.numeric(unlist(elt2[lag,]))
            row <- as.integer(hr+lag-1)
            set(air.mass,row,1:nc,mc.air+k1a*t1+k2a*t2)
            set(sur.mass,row,1:nc,mc.sur+k1s*t1+k2s*t2)
          }  
        }  
      }
    } else {
      air.mass  <- as.data.table(matrix(rep(mc.air,each=8736),nrow=8736))
      sur.mass  <- as.data.table(matrix(rep(mc.air,each=8736),nrow=8736))
    }  
    day.air  <- colMeans(array(as.matrix(air.mass),c(24,364,nc)),dim=1)
    day.sur  <- colMeans(array(as.matrix(sur.mass),c(24,364,nc)),dim=1)
    day.win  <- day.air*flows[1]$aer.out
    day.was  <- day.air*flows[1]$cln.air + day.sur*flows[1]$cln.sur
    fug.day  <- as.data.table(cbind(1:364,rep(0,364),day.air,day.sur,day.win,day.was))
    setnames(fug.day,c("daynum","hour",str_c(c(rep("air",nc),rep("sur",nc),rep("win",nc),rep("was",nc)),rep(1:nc,4))))
    padding  <- matrix(0,nrow=8736,ncol=2*nc)
    fug.hour <- as.data.table(cbind(unlist(lapply(1:364,rep,24)),rep(1:24,364),air.mass,sur.mass,padding))
    setnames(fug.hour,c("daynum","hour",str_c(c(rep("air",nc),rep("sur",nc),rep("win",nc),rep("was",nc)),rep(1:nc,4))))
    fugs     <- rbind(fug.day,fug.hour)
    if (g$prog=="y") cat("\n  Evaluating fugacity concentrations complete")
    return(fugs)
  }
  
  
  
  
  # eval.fug.concs.st solves the fugacity equations using a finite time step 
  
  eval.fug.concs.st = function(chem.release,flows,days,steps) {
    nc          <- nrow(flows)
    air.mass    <- as.data.table(matrix(0,nrow=8736,ncol=nc))
    sur.mass    <- as.data.table(matrix(0,nrow=8736,ncol=nc))
    waste.mass  <- as.data.table(matrix(0,nrow=364,ncol=nc))
    window.mass <- as.data.table(matrix(0,nrow=364,ncol=nc))
    yair        <- select_vars(names(chem.release),starts_with("air"))
    ysur        <- select_vars(names(chem.release),starts_with("sur"))
    df          <- as.data.frame(chem.release)
    new.air     <- as.data.table(df[yair])     # chemical mass in mg 
    new.sur     <- as.data.table(df[ysur])     # chemical mass in mg 
    air <- rep(0,nc) 
    sur <- rep(0,nc)
    was <- rep(0,nc)
    win <- rep(0,nc)
    fac   <- 24*steps
    fromair <- 1-exp(-flows$a/fac) 
    toair   <- 1-exp(-flows$b/fac)
    tosur   <- 1-exp(-flows$c/fac)
    fromsur <- 1-exp(-flows$d/fac)
    towin   <- 1-exp(-flows$aer.out/fac)
    cl.air  <- 1-exp(-flows$cln.air/fac)
    cl.sur  <- 1-exp(-flows$cln.sur/fac)
    for (day in 1:days) {
      win  <- win*0
      was  <- was*0
      doff <- 24*(day-1)
      for (hr in 1:24) {
        i <- as.integer(doff+hr)
        air <- air + new.air[i]
        sur <- sur + new.sur[i]
        for (j in 1:steps) {
          air <- air - air*fromair + sur*toair
          sur <- sur - sur*fromsur + air*tosur
          win <- win + air*towin
          was <- was + air*cl.air + sur*cl.sur
        }
        set(air.mass,i,1:nc,air)
        set(sur.mass,i,1:nc,sur)
      }  
      set(window.mass,day,1:nc,win)
      set(waste.mass,day,1:nc,was)
    }
    return(list(air.mass,sur.mass.waste.mass,window.mass))
  }
  
  
  
  # eval.hm rate sets the hand-to-mouth rate constant for the primary person
  # based on results from SHEDS soil-dust model runs (up to age 20) 
  # top age group extended to all adults
  
  eval.hm.rate = function(prime,q) {
    age <- prime$age
    if (age==0) {
      gm  <- 0.20
      gsd <- 1.8
    } else if (age==1) {
      gm  <- 0.23
      gsd <- 1.7
    } else if (age==2) {
      gm  <- 0.19
      gsd <- 1.9
    } else if (age>=3 & age<=5) {
      gm  <- 0.15
      gsd <- 1.9
    } else if (age>=6 & age<=10) {
      gm  <- 0.09
      gsd <- 2.0
    } else if (age>=11 & age<=15) {
      gm  <- 0.05
      gsd <- 2.4
    } else if (age>=16) {
      gm  <- 0.02
      gsd <- 3.1
    } 
    gm <- gm/2    # convert from fraction to rate, assuming 2 hrs average before hand washing
    kh <- distrib("logn",gm,gsd,q=q)
    return(kh)
  }
  
  
  
  # eval.house.props samples distributions for house characteristics relating to fugacity calculations
  
  eval.house.props <- function(fug.hvars,person.data,q.house) {
    if (g$prog=="y") cat("\n  Evaluating house props...")
    area.sur     <- person.data$unitsf/10.7639    # convert from ft2 to m2
    f            <- fug.hvars[fug.hvars$varname=="aer.out"]
    aer.out      <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$aer.out)
    f            <- fug.hvars[fug.hvars$varname=="height"]
    height       <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$height)
    f            <- fug.hvars[fug.hvars$varname=="lg.carb.f"]
    lg.carb.f    <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$lg.carb.f)
    f            <- fug.hvars[fug.hvars$varname=="lg.clean.air"]
    lg.clean.air <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$lg.clean.air)
    f            <- fug.hvars[fug.hvars$varname=="lg.clean.sur"]
    lg.clean.sur <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$lg.clean.sur)
    f            <- fug.hvars[fug.hvars$varname=="lg.depos"]
    lg.depos     <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$lg.depos)
    f            <- fug.hvars[fug.hvars$varname=="lg.load.air"]
    lg.load.air  <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$lg.load.air)
    f            <- fug.hvars[fug.hvars$varname=="lg.load.sur"]
    lg.load.sur  <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$lg.load.sur)
    lg.load.sur  <- 1E4*lg.load.sur     # convert from ug/cm2 to ug/m2
    f            <- fug.hvars[fug.hvars$varname=="lg.resus"]
    lg.resus     <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$lg.resus)
    f            <- fug.hvars[fug.hvars$varname=="sm.carb.f"]
    sm.carb.f    <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$sm.carb.f)
    f            <- fug.hvars[fug.hvars$varname=="sm.clean.air"]
    sm.clean.air <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$sm.clean.air)
    f            <- fug.hvars[fug.hvars$varname=="sm.clean.sur"]
    sm.clean.sur <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$sm.clean.sur)
    f            <- fug.hvars[fug.hvars$varname=="sm.depos"]
    sm.depos     <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$sm.depos)
    f            <- fug.hvars[fug.hvars$varname=="sm.load.air"]
    sm.load.air  <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$sm.load.air)
    f            <- fug.hvars[fug.hvars$varname=="sm.load.sur"]
    sm.load.sur  <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$sm.load.sur)
    sm.load.sur  <- 1E4*sm.load.sur    # convert from ug/cm2 to ug/m2
    f            <- fug.hvars[fug.hvars$varname=="sm.resus"]
    sm.resus     <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$sm.resus)
    f            <- fug.hvars[fug.hvars$varname=="temp"]
    temp         <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$temp)
    f            <- fug.hvars[fug.hvars$varname=="thick.bou"]
    thick.bou    <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$thick.bou)
    f            <- fug.hvars[fug.hvars$varname=="thick.sur"]
    thick.sur    <- distrib(f$form,f$par1,f$par2,f$par3,f$par4,f$lower.trun,f$upper.trun,f$resamp,q=q.house$thick.sur)
    house        <- person.data$house
    house.props  <- as.data.table(cbind(area.sur,aer.out,height,lg.carb.f,lg.clean.air,lg.clean.sur,
                                        lg.depos,lg.load.air,lg.load.sur,lg.resus,sm.carb.f,sm.clean.air,sm.clean.sur,
                                        sm.depos,sm.load.air,sm.load.sur,sm.resus,temp,thick.bou,thick.sur,house))
    setnames(house.props,c("area.sur","aer.out","height","lg.carb.f","lg.clean.air","lg.clean.sur",
                           "lg.depos","lg.load.air","lg.load.sur","lg.resus","sm.carb.f","sm.clean.air","sm.clean.sur",
                           "sm.depos","sm.load.air","sm.load.sur","sm.resus","temp","thick.bou","thick.sur","house"))
    return(house.props)
  }
  
  
  
  
  # eval.indirect calculates the indirect exposures on a daily basis
  
  eval.indirect = function(d,fug.hour,dermal.rates,hp,nc,prime,fug.cvars,chem.list) {
    # first, find the status for the primary person at each hour
    # 3 options:  sleep, awake, out
    # each hour of yearF has # minutes in each: min.sleep, min.awake, min.out
    # start.min and end.min refer to status, not to product use
    df <- as.data.frame(d)
    df$dayofweek  <- 1 + (df$daynum-1) %% 7
    df$weekend    <- 0
    df$weekend[df$dayofweek>5] <- 1
    df$start.min  <- 60*(df$hournum-1)+round(60*(1+df$start-df$hour))
    df$end.min    <- lead(df$start.min)-1
    df$end.min[is.na(df$end.min)]<- 8736*60
    # activites: -1=idle, -2=PUC, 1,2=commute, 3,4,5=eat B,D,Lunch, 6=sleep, 7=work
    # status   : 1=sleep, 2=awake, 3=away from home
    df$status <- 2
    df$status[df$activity_code==6] <- 1 
    df$status[df$activity_code==7] <- 3
    df$status[df$activity_code==1] <- 3
    df$status[df$activity_code==2] <- 3
    df$status[df$activity_code==5] <- 3           # assume lunch is away from home, other meals at home
    status.min <- rep(1,8736*60)                  # status.min has value for each minute of year
    for (i in 1:nrow(df)) {
      status.min[df$start.min[i]:df$end.min[i]] <- df$status[i]
    }
    
    fh <- as.data.frame(fug.hour)
    fh$min.sleep <- rep(0,8736)
    fh$min.awake <- rep(0,8736)
    fh$min.out   <- rep(0,8736)
    dt           <- as.data.table(fh)
    rates        <- as.data.table(dermal.rates)[dermal.rates$source.id=="Indirect"]
    x <- matrix(status.min,nrow=60,ncol=8736)
    for (i in 1:8736){
      y <- plyr::count(x[,i])
      sleep <- y$x==1
      awake <- y$x==2
      out   <- y$x==3
      if(any(sleep)) set(dt,i,"min.sleep",y$freq[sleep])
      if(any(awake)) set(dt,i,"min.awake",y$freq[awake])
      if(any(out))   set(dt,i,"min.out",  y$freq[out])
    }
    
    fug <- as.data.frame(dt)
    air.max <- matrix(0,nrow=364,ncol=nc)
    air.exp <- matrix(0,nrow=364,ncol=nc)
    air.mass<- matrix(0,nrow=364,ncol=nc)
    air.abs <- matrix(0,nrow=364,ncol=nc)
    der.exp <- matrix(0,nrow=364,ncol=nc)
    der.max <- matrix(0,nrow=364,ncol=nc)
    der.abs <- matrix(0,nrow=364,ncol=nc)
    ing.exp <- matrix(0,nrow=364,ncol=nc)
    ing.abs <- matrix(0,nrow=364,ncol=nc)
    tc.hands <- 300                           # transfer coefficient for hands is 300 cm2/hr (assumption)
    hand.area <-  0.05*prime$skin.area        # assumed hand skin area is 5% of total 
    vol.cloud <- 2                            # personal cloud has volume of 2 m3
    aer.cloud <- 10                           # personal cloud exchanges its air 10 times per hour.
    
    for (c in 1:nc) {
      survar       <- str_c("sur",c)
      airvar       <- str_c("air",c)
      if (sum(fug[survar],fug[airvar])>0) {
        ka <- max(0.00001,rates$ka[rates$dtxsid==chem.list[c]])
        ks <- max(0.00001,rates$ks[rates$dtxsid==chem.list[c]])
        kh <- max(0.00001,rates$kh[rates$dtxsid==chem.list[c]])
        kr <- max(0.00001,rates$kr[rates$dtxsid==chem.list[c]])
        hw <- max(0.00001,rates$hw[rates$dtxsid==chem.list[c]])
        dur.indirect <- 8/hw                                              # indirect exposure lasts (8 hr /handwashes per day)
        f.hands.lost <- 1-exp(-(ka+ks+kh+kr)*dur.indirect)                # amount lost from hands
        sur.eff      <- as.vector(unlist(fug[survar]))
        sur.conc     <- as.vector(sur.eff/(10000*hp$area.sur))            # units are converted to [mg/cm2]
        f.avail      <- 0.5                                               # assumed availability fraction
        hands.exp    <- sur.conc * tc.hands * f.avail * fug$min.awake/60  # units are [mg] for current hour
        hands.ldg    <- hands.exp/hand.area
        hands.abs    <- hands.exp * ks/(ka+ks+kh+kr) * f.hands.lost       # assumes effect is all on current hour
        hands.air    <- hands.exp * ka/(ka+ks+kh+kr) * f.hands.lost
        hands.ing    <- hands.exp * kh/(ka+ks+kh+kr) * f.hands.lost
        hands.aconc  <- hands.air / (vol.cloud*aer.cloud*dur.indirect)    # personal cloud conc. from hands
        
        air.eff      <- as.vector(unlist(fug[airvar]))
        air.conc     <- air.eff/(hp$area.sur*hp$height) + hands.aconc     # add house term and personal cloud
        inhal.exp    <- air.conc * (fug$min.sleep+fug$min.awake)/60
        inhal.mass1  <- air.conc * fug$min.sleep/60 * prime$basal.vent/24         # convert m3/day to m3/hr
        inhal.mass2  <- air.conc * fug$min.awake/60 * prime$basal.vent/24 * 2.2   # assume MET=2.2 when awake
        inhal.mass   <- inhal.mass1 + inhal.mass2
        inhal.abs    <- inhal.mass * 0.16                                 # assume 16% lung absorption
        inges.abs    <- hands.ing * fug.cvars$fabs[c]
        m.h.exp      <- matrix(hands.exp,nrow=24,ncol=364)
        m.h.ldg      <- matrix(hands.ldg,nrow=24,ncol=364)
        m.h.abs      <- matrix(hands.abs,nrow=24,ncol=364)
        m.h.ing      <- matrix(hands.ing,nrow=24,ncol=364)
        m.a.exp      <- matrix(inhal.exp,nrow=24,ncol=364)
        m.a.conc     <- matrix(air.conc,nrow=24,ncol=364)
        m.a.mass     <- matrix(inhal.mass,nrow=24,ncol=364) 
        m.a.abs      <- matrix(inhal.abs,nrow=24,ncol=364)
        m.g.abs      <- matrix(inges.abs,nrow=24,ncol=364)
        for (i in 1:364) {
          der.exp[i,c] <- sum(m.h.exp[,i])
          der.max[i,c] <- max(m.h.ldg[,i])
          der.abs[i,c] <- sum(m.h.abs[,i])
          ing.exp[i,c] <- sum(m.h.ing[,i])
          air.exp[i,c] <- sum(m.a.exp[,i])/24
          air.max[i,c] <- max(m.a.conc[,i])
          air.mass[i,c]<- sum(m.a.mass[,i])
          air.abs[i,c] <- sum(m.a.abs[,i])
          ing.abs[i,c] <- sum(m.g.abs[,i])
        }   # end day loop
      }     # end if c present
    }       # end c loop  
    nam.de    <- str_c(rep("ind.derm.exp",nc),1:nc)
    nam.dm    <- str_c(rep("ind.derm.max",nc),1:nc)
    nam.da    <- str_c(rep("ind.derm.abs",nc),1:nc)
    nam.ihe   <- str_c(rep("ind.inhal.exp",nc),1:nc)
    nam.ihm   <- str_c(rep("ind.inhal.max",nc),1:nc)
    nam.ing   <- str_c(rep("ind.ingest.exp",nc),1:nc)
    nam.am    <- str_c(rep("ind.inhal.mass",nc),1:nc)
    nam.aa    <- str_c(rep("ind.inhal.abs",nc),1:nc)
    nam.ga    <- str_c(rep("ind.ingest.abs",nc),1:nc)
    indirect <- data.table(1:364,der.exp,der.max,der.abs,air.exp,air.max,air.mass,air.abs,ing.exp,ing.abs)
    setnames(indirect,c("daynum",nam.de,nam.dm,nam.da,nam.ihe,nam.ihm,nam.am,nam.aa,nam.ing,nam.ga))
    if (g$prog=="y") cat("\n  Evaluating indirect exposures complete")
    return(indirect)
  }
  
  
  
  
  # eval.lcia calculates the PIF fractions for each household
  
  eval.lcia = function(puc.list,chem.list,prod.chem,diary,annual,house.num,nc) {
    lcia <- NULL
    if(nrow(annual)>0) {
      lcia <- data.table(rep(puc.list,nc),chem.list,matrix(0,nrow=nc,ncol=11))
      setnames(lcia,c("source.id","dtxsid","prod.mass","frac.chem","chem.mass","pif.derm","pif.inhal","pif.ingest","mass.air",
                      "mass.water","mass.land","age","household"))
      num.persons <- length(unique(diary$p))
      for (c in 1:length(chem.list)) {
        mass      <- sum(as.numeric(diary$mass[diary$source.id==puc.list & diary$primary==1]))
        totmass   <- sum(as.numeric(diary$mass[diary$source.id==puc.list]))
        year      <- annual[annual$dtxsid==chem.list[c]]
        lcia[c,3] <- mass
        lcia[c,4] <- select(prod.chem,c+3)
        chem.mass <- mass * select(prod.chem,c+3)
        lcia[c,5] <- chem.mass
        factor    <- num.persons*mass/totmass
        if(chem.mass>0) { 
          lcia[c,6]  <- 364*(year$dir.derm.exp   + factor * year$ind.derm.exp)  / (1000*chem.mass)
          lcia[c,7]  <- 364*(year$dir.inhal.mass + factor * year$ind.inhal.mass)/ (1000*chem.mass)
          lcia[c,8]  <- 364*(year$dir.ingest.exp + factor * year$ind.ingest.exp)/ (1000*chem.mass)
          lcia[c,9]  <- 364*year$out.air * mass/totmass / 1000                         # convert back to grams
          lcia[c,10] <- 364*year$drain * mass/totmass  / 1000                          # convert back to grams
          lcia[c,11] <- 364*(year$out.sur+year$waste) * mass/totmass / 1000            # convert back to grams
        }
        lcia[c,12] <- unique(diary$age[diary$primary==1])
        lcia[c,13] <- house.num
      }  
    }  
    return(lcia)
  }
  
  
  
  # eval.lcia.avg averages the LCIA values across households
  
  eval.lcia.avg = function(all.lcia,puc.list,chem.list,nc) {
    lcia.avg <- data.table(rep(puc.list,nc),chem.list,matrix(0,nrow=nc,ncol=13))
    setnames(lcia.avg,c("source.id","dtxsid","mass.frac.chem.puc","mass.puc.use.adult","pif.derm.adult","pif.inhal.adult","pif.ingest.adult",
                        "mass.puc.use.child","pif.derm.child","pif.inhal.child","pif.ingest.child",
                        "chem.mass","mass.tot.air","mass.tot.water","mass.tot.land"))
    for (c in 1:length(chem.list)) {
      lcia.chem      <- all.lcia[all.lcia$dtxsid==chem.list[c]]
      lcia.avg[c,3]  <- mean(lcia.chem$frac.chem)
      lcia.avg[c,4]  <- mean(lcia.chem$prod.mass[lcia.chem$age>=18])
      lcia.avg[c,5]  <- mean(lcia.chem$pif.derm[lcia.chem$age>=18])
      lcia.avg[c,6]  <- mean(lcia.chem$pif.inhal[lcia.chem$age>=18])
      lcia.avg[c,7]  <- mean(lcia.chem$pif.ingest[lcia.chem$age>=18])
      lcia.avg[c,8]  <- mean(lcia.chem$prod.mass[lcia.chem$age<18])
      lcia.avg[c,9]  <- mean(lcia.chem$pif.derm[lcia.chem$age<18])
      lcia.avg[c,10] <- mean(lcia.chem$pif.inhal[lcia.chem$age<18])
      lcia.avg[c,11] <- mean(lcia.chem$pif.ingest[lcia.chem$age<18])
      lcia.avg[c,12] <- mean(lcia.chem$chem.mass)
      lcia.avg[c,13] <- mean(lcia.chem$mass.air)
      lcia.avg[c,14] <- mean(lcia.chem$mass.water)
      lcia.avg[c,15] <- mean(lcia.chem$mass.land)
    }
    return(lcia.avg)
  }
  
  
  
  #eval.prime detemines certain personal variables
  
  eval.prime = function(persons,ventilation,q) {
    prime            <- persons[persons$primary==1]
    prime$hand.wash  <- distrib("logn",3.74,2.63,lt=1,ut=12,q=as.numeric(q$hand.wash))
    prime$hand.mouth <- eval.hm.rate(prime,as.numeric(q$hand.mouth))
    if (!exists("basal.vent",prime)) {
      b     <- ventilation[minage<=prime$age & prime$age<=maxage & prime$sex==sex]
      qtrun <- 0.02+0.96*q$basal.vent                               # truncate 2% at each end
      prime$basal.vent <- distrib("norm",b$mean_ve,b$std_ve,q=qtrun)
    }
    return(prime)
  }
  
  
  
  # eval.prod.chem produces a table of chemical masses and compartment fractions for each PUC
  
  eval.prod.chem = function(pucs,puc.list,brand.list,chem.list,chem.fracs,release.fracs,q,house.num) {
    y  <- pucs$source.id
    nc <- length(chem.list)
    f  <- matrix(0,nrow=length(y),ncol=nc)
    prod <- vector("integer",length(y))
    form <- vector("integer",length(y))
    if (length(y)>0) {
      for (i in 1:length(y)) {
        p <- y[i]
        qname1 <- str_c("prod.id",which(puc.list==p))
        qname2 <- str_c("form.id",which(puc.list==p))
        z <- chem.fracs[chem.fracs$source.id==p]
        nprod   <- brand.list$brands[brand.list$puc==p][[1]]
        prod[i] <- nprod[ceiling(q[[qname1]]*length(nprod))]
        nform   <- unique(z$formulation_id[z$product_id==prod[i]])
        choice <- NULL
        if(length(nform)>0) {
          form[i] <- nform[ceiling(q[[qname2]]*length(nform))]
          choice  <- z[z$product_id==prod[i] & z$formulation_id==form[i]]
        }  
        for (j in 1:length(chem.list)) {
          cc <- choice[choice$dtxsid==chem.list[j]]
          if (!is.null(cc)) {if (nrow(cc)==1) f[i,j] <- cc$weight_fraction  }
        }
      }
    }
    t <- as.data.table(data.frame(y,prod,form,f))
    t2<-t
    setnames(t2,c("source.id","product_id","formulation_id",chem.list))
    t2long<-melt(t2, c("source.id","product_id","formulation_id"))
    prod_chem_ids <-t2long[!t2long$value ==0]
    setnames(prod_chem_ids,c("puc","product_id","formulation_id","chemical","mass.fraction"))
    write.csv(prod_chem_ids,paste0("output/S2D/",run.name,"/Prod_chem","/Prod_chem_",house.num,".csv"),row.names=FALSE) 
    
    setnames(t,c("source.id","product_id","formulation_id",str_c(rep("X",nc),rep(1:nc))))
    prod.chem <- left_join(t,release.fracs,by="source.id")
    prod.chem$fia[is.na(prod.chem$fia)] <- 0
    prod.chem$fis[is.na(prod.chem$fis)] <- 0
    prod.chem$foa[is.na(prod.chem$foa)] <- 0
    prod.chem$fos[is.na(prod.chem$fos)] <- 0
    prod.chem$fsk[is.na(prod.chem$fsk)] <- 0
    prod.chem$fws[is.na(prod.chem$fws)] <- 0
    prod.chem$fdr[is.na(prod.chem$fdr)] <- 0
    prod.chem$fap[is.na(prod.chem$fap)] <- 0
    prod.chem$fsl[is.na(prod.chem$fsl)] <- 0
    prod.chem$fgi[is.na(prod.chem$fgi)] <- 0
    prod.chem$fhn[is.na(prod.chem$fhn)] <- 0
    if (g$prog=="y") cat("\n  Evaluating product chemicals complete")
    return(prod.chem)
  }
  
  
  
  
  # eval.pucs.areas determines the affected skina areas for each puc
  
  eval.pucs.areas = function(pucs,skin.areas,q) {
    x <- left_join(pucs,select(skin.areas,source.id,adult,child,hand_variability,body_variability),by="source.id")
    mode(x$mass) <- "numeric"
    x$f.adult <- 0.01
    x$f.child <- 0.01
    if (nrow(x)>0) {
      x$adult[is.na(x$adult)] <- 0.01
      x$child[is.na(x$child)] <- 0.01
      x$hand_variability[is.na(x$hand_variability)] <- "low"
      x$body_variability[is.na(x$body_variability)] <- "low"
      x$hand_mult[x$hand_variability=="low"]    <- 12
      x$hand_mult[x$hand_variability=="medium"] <- 4
      x$hand_mult[x$hand_variability=="high"]   <- 1
      x$body_mult[x$body_variability=="low"]    <- 12
      x$body_mult[x$body_variability=="medium"] <- 4
      x$body_mult[x$body_variability=="high"]   <- 1
      for (i in 1:nrow(x)) {
        qname <- str_c("skin.frac",i)
        if (x$hand_variability[i]!="none") {
          cadult <- x$hand_mult[i]
          dadult <- cadult/x$adult[i]-cadult
          cchild <- x$hand_mult[i]
          dchild <- cchild/x$child[i]-cchild
          x$f.adult[i] <- qbeta(q[[qname]],cadult,dadult)*0.05
          x$f.child[i] <- qbeta(q[[qname]],cchild,dchild)*0.05
        }
        if (x$body_variability[i]!="none") {
          cadult <- x$body_mult[i]
          dadult <- cadult/x$adult[i]-cadult
          cchild <- x$body_mult[i]
          dchild <- cchild/x$child[i]-cchild
          x$f.adult[i] <- x$f.adult[i] + qbeta(q[[qname]],cadult,dadult)*0.95
          x$f.child[i] <- x$f.child[i] + qbeta(q[[qname]],cchild,dchild)*0.95
        }
      }  
    }
    x$h                <- NA
    y                  <- x[1]
    y$source.id        <- "Indirect"
    y$product_type     <- "none"
    y$indoor_outdoor   <- "I"
    y$hand.dur         <- 0
    y$mass             <- 0
    y$code             <- "IND"
    y$met              <- 2.2
    y$spray            <- FALSE
    y$adult            <- 0.04
    y$child            <- 0.04
    y$hand_variability <- "none"
    y$body_variability <- "none"
    y$f.adult          <- 0.04
    y$f.child          <- 0.04
    y$h                <- 1E-6 + 9E-6*q$thick.indir
    x                  <- rbind(x,y)
    return(x)
  }
  
  
  
  
  # eval.puc.wipe.rinse merges the puc information with the removal fractions.
  
  eval.puc.wipe.rinse = function(pucs,removal.fracs) {
    x <- as.data.frame(removal.fracs)
    y <- data.table(x$source.id,x$fhands,x$fbody,x$frinseh,x$frinseb,x$fwipeh,x$fwipeb,x$frinses,x$fwipes)
    setnames(y,c("source.id","fhands","fbody","frinseh","frinseb","fwipeh","fwipeb","frinses","fwipes"))
    z <- y[y$source.id %in% pucs$source.id]
    setkey(z,source.id)
    setkey(pucs,source.id)
    puc.wipe.rinse <- join(pucs,z,by="source.id")
    puc.wipe.rinse[is.na(puc.wipe.rinse)] <- 0.5
    return(puc.wipe.rinse)
  }
  
  
  
  # eval.release.fracs lists the compartmental fractional releases for each PUC
  
  eval.release.fracs = function(puc.types,compart.fracs) {
    if (g$prog=="y") cat("\n  Evaluating release fractions...")
    x <- as.data.table(left_join(puc.types,compart.fracs,by="code"))
    setorder(x,source.id)
    release.fracs <- x
    return(release.fracs)
  }
  
  
  
  # eval.release.times generates a calendar of chemical releases
  
  eval.release.times = function(chem.release) {
    release.hrs  <- as.vector(1:8736)[chem.release$rel.tot>1E-9]    # at least 10^-9 milligrams released
    indoor.hrs   <- as.vector(1:8736)[chem.release$rel.ind>1E-9]
    nr           <- length(release.hrs)
    xr           <- c(release.hrs,8737)
    release.gaps <- xr[2:(nr+1)]-xr[1:nr]-1
    ni           <- length(indoor.hrs)
    xi           <- c(indoor.hrs,8737)
    indoor.gaps  <- xi[2:(ni+1)]-xi[1:ni]-1
    if (g$prog=="y") cat("\n  Evaluating release times complete")
    return(list(release.hrs,release.gaps,indoor.hrs,indoor.gaps))
  }
  
  
  
  # eval.summary writes daily and annual output files for one household
  
  eval.summary = function(direct,indirect,env.impact,run.name,house.num,nc,chem.list,chem.totals) {
    var.list  <- c("dir.derm.exp","dir.derm.max","dir.derm.abs","dir.inhal.exp","dir.inhal.max","dir.inhal.mass",
                   "dir.inhal.abs","dir.ingest.exp","dir.ingest.abs","dir.release","ind.derm.exp","ind.derm.max",
                   "ind.derm.abs","ind.inhal.exp","ind.inhal.max","ind.inhal.mass","ind.inhal.abs",
                   "ind.ingest.exp","ind.ingest.abs","out.sur","out.air","drain","waste")
    ndays <- 364
    nvars <- length(var.list)
    if(!is.null(nrow(direct)))   { x <- select(direct,-daynum) } else {
      x <- data.table(matrix(0,nrow=ndays,ncol=10*nc))
      dir.names <- select_vars(var.list,starts_with("dir"))
      setnames(x, unlist(lapply(dir.names,str_c,1:nc)))
    }
    if(!is.null(nrow(indirect)))   { y <-select(indirect,-daynum) } else {
      y <- data.table(matrix(0,nrow=ndays,ncol=9*nc))
      ind.names <- select_vars(var.list,starts_with("ind"))
      setnames(y, unlist(lapply(ind.names,str_c,1:nc)))
    }
    if(!is.null(nrow(env.impact)) & ncol(env.impact)>1) { z <-select(env.impact,-daynum) } else {
      z <- data.table(matrix(0,nrow=ndays,ncol=4*nc))
      env.names <- var.list[(nvars-3):nvars]
      setnames(z, unlist(lapply(env.names,str_c,1:nc)))
    }  
    daily.base <- data.table(x,y,z) 
    daily.data <- matrix(0,nrow=ndays*nc,ncol=nvars)
    daily.info <- data.table(rep(house.num,ndays*nc),rep(chem.list,each=ndays),rep(1:ndays,nc))
    setnames(daily.info,c("household","dtxsid","daynum"))
    for (c in 1:nc) {
      rows <- daily.info$dtxsid==chem.list[c]
      for (j in 1:nvars) {
        z <- daily.base[[paste0(var.list[j],c)]]
        if(!is.null(z)) daily.data[rows,j] <- z
      }
    }
    daily <- data.table(daily.info,daily.data)
    setnames(daily,c(names(daily.info),var.list))
    write.csv(daily,paste0("output/S2D/",run.name,"/Daily","/Daily_",house.num,".csv"),row.names=FALSE)
    # now do the annual summary
    annual.base <- colMeans(daily.base)
    annual.data <- matrix(0,nrow=nc,ncol=nvars)
    annual.info <- data.table(rep(house.num,nc),chem.list)
    setnames(annual.info,c("household","dtxsid"))
    for (c in 1:nc){
      for (j in 1:nvars) {
        z <- annual.base[[paste0(var.list[j],c)]]
        if(!is.null(z)) annual.data[c,j] <- z
      }
    }
    annual.data <- as.data.table(annual.data)
    setnames(annual.data,var.list)
    totals <- data.table(matrix(unlist(chem.totals),nrow=nc,ncol=1))
    setnames(totals,"total.used")
    annual <- data.table(annual.info,totals,annual.data)
    write.csv(annual,paste0(g$out,"/Annual/House_",house.num,".csv"),row.names=FALSE)
    return(annual)
  }
  
  
  
  # eval.use.chem gives the chemical masses released in total by each product use
  
  eval.use.chem = function(pucs,prod.chem,nc) {
    pm       <- pucs$mass*1000                     # use.chem has units of (mg per use)
    y1       <- select_vars(names(prod.chem),starts_with("X"))
    z1       <- as.matrix(as.data.frame(prod.chem)[y1])
    use.chem <- matrix(0,nrow=nrow(pucs),ncol=nc)
    for (prod in 1:nrow(pucs)) {
      use.chem[prod,]    <- z1[prod,]*pm[prod]
    }    
    use.chem <- as.data.table(use.chem)
    setnames(use.chem,str_c(rep("tot",nc),1:nc))
    if (g$prog=="y") cat("\n  Evaluating use.chem complete")
    return(use.chem)
  }
  
  
  
  # eval.use.data creates a matrix of chemical masses in 9 compartments, for each PUC
  eval.use.data = function(pucs,prod.chem,fug.cvars,dermal.rates,puc.wipe.rinse,nc,house.num) {
    pm  <- pucs$mass*1000                                                               # chemical mass from grams to milligrams [mg]
    df  <- as.data.frame(prod.chem)                                               
    y1  <- select_vars(names(prod.chem),starts_with("X"))                               # mass fractions by chemical
    y2  <- select_vars(names(prod.chem),starts_with("f"),exclude="formulation_id")      # mass fractions by compartment
    z1  <- as.matrix(as.data.frame(prod.chem)[y1])
    z2  <- as.matrix(as.data.frame(prod.chem)[y2])
    use.data <- array(0,dim=c(nrow(pucs),3,length(y2),nc))
    for (prod in 1:nrow(pucs)) {
      use.data[prod,1,,] <- (z2[prod,] %o% z1[prod,])*pm[prod]                          # mass of chemicals in each of 11 compartments initially
      vp     <- fug.cvars$vapor*760/101325                                              # convert from Pascals to Torr for use in ke formula
      ke     <- (log(10)/145/60)*(fug.cvars$molwt*vp)^0.9546                            # ke in [1/min] is volatilization rate from indoor surface
      rates  <- dermal.rates[source.id==pucs$source.id[prod]] 
      ka     <- as.data.frame(rates)$ka/60                                       # Ka in [1/min] is volatilization rate from skin
      f.sur  <- 1-exp(-ke*pucs$hand.dur[prod])                                          # fraction of surface chemical lost to air
      f.skin <- 1-exp(-ka*pucs$hand.dur[prod])                                          # fraction of chemical on skin lost to air
      toair1 <- use.data[prod,1,y2=="fis",]*f.sur                                       # hand.dur in [min] - surface to air during indoor use
      toair2 <- use.data[prod,1,y2=="fos",]*f.sur                                       # hand.dur in [min] - surface to air during outdoor use
      toair3 <- use.data[prod,1,y2=="fsk",]*f.skin                                      # loss from skin to air during use
      if (pucs[prod]$spray==TRUE) {                                         
        tosur1 <- use.data[prod,1,y2=="fia",]*0.6                                      # 60% of spray settles during use  
        tosur2 <- use.data[prod,1,y2=="foa",]*0.6                                      # 60% of spray settles during use
        use.data[prod,1,y2=="fia",] <- use.data[prod,1,y2=="fia",] - tosur1            # adjust tier 1 for settling indoors
        use.data[prod,1,y2=="fis",] <- use.data[prod,1,y2=="fis",] + tosur1            # adjust tier 1 for settling indoors
        use.data[prod,1,y2=="foa",] <- use.data[prod,1,y2=="foa",] - tosur2            # adjust tier 1 for settling outdoors
        use.data[prod,1,y2=="fos",] <- use.data[prod,1,y2=="fos",] + tosur2            # adjust tier 1 for settling outdoors
      }    
      # Now start tier 2 calculations
      use.data[prod,2,,] <- use.data[prod,1,,]                                          # tier 2 = tier 1
      if(prod.chem$ind[prod]==1) {                                                      # indoor transfers    
        use.data[prod,2,y2=="fia",] <- use.data[prod,2,y2=="fia",] + toair1             # surface to air transfer
        use.data[prod,2,y2=="fia",] <- use.data[prod,2,y2=="fia",] + toair3             # skin to air transfer 
        use.data[prod,2,y2=="fis",] <- use.data[prod,2,y2=="fis",] - toair1             # loss from surface
        use.data[prod,2,y2=="fsk",] <- use.data[prod,2,y2=="fsk",] - toair3             # loss from skin
      }
      if(prod.chem$ind[prod]==0) {                                                      # outdoor transfers    
        use.data[prod,2,y2=="foa",] <- use.data[prod,2,y2=="foa",] + toair2             # surface to air transfer
        use.data[prod,2,y2=="foa",] <- use.data[prod,2,y2=="foa",] + toair3             # skin to air transfer 
        use.data[prod,2,y2=="fos",] <- use.data[prod,2,y2=="fos",] - toair2             # loss from surface
        use.data[prod,2,y2=="fsk",] <- use.data[prod,2,y2=="fsk",] - toair3             # loss from skin
      }
      # Tier 2 reflects end of use, before rinse off and wipe off
      # Tier 3 reflects the post-use values after rinse off and wipe off
      use.data[prod,3,,] <- use.data[prod,2,,]                                          # tier 3 starts equal to tier 2
      if (pucs[prod]$spray==TRUE) {              
        tosur3 <- pmin(use.data[prod,1,y2=="fia",],use.data[prod,3,y2=="fia",])        # remainder of spray settles when use ends  
        tosur4 <- pmin(use.data[prod,1,y2=="foa",],use.data[prod,3,y2=="foa",])        # remainder of spray settles when use ends
        use.data[prod,3,y2=="fia",] <- use.data[prod,3,y2=="fia",] - tosur3            # adjust tier 3 for settling indoors
        use.data[prod,3,y2=="fis",] <- use.data[prod,3,y2=="fis",] + tosur3            # adjust tier 3 for settling indoors
        use.data[prod,3,y2=="foa",] <- use.data[prod,3,y2=="foa",] - tosur4            # adjust tier 3 for settling outdoors
        use.data[prod,3,y2=="fos",] <- use.data[prod,3,y2=="fos",] + tosur4            # adjust tier 3 for settling outdoors
      }    
      removals <- puc.wipe.rinse[puc.wipe.rinse$source.id==pucs$source.id[prod]]        # puc-specific removal terms
      todrains <- use.data[prod,2,y2=="fis",]*removals$frinses                          # loss from surfaces due to rinse off
      hands1   <- use.data[prod,2,y2=="fsk",]*removals$fhands                           # skin loading on hands before rinse/wipe
      body1    <- use.data[prod,2,y2=="fsk",]*removals$fbody                            # skin loading on body before rinse/wipe
      todrainh <- hands1*removals$frinseh                                               # amount rinsed off hands
      todrainb <- body1*removals$frinseb                                                # amount rinsed off body
      towastes <- use.data[prod,2,y2=="fis",]*removals$fwipes                           # loss from surfaces due to wipe off 
      towasteh <- hands1*removals$fwipeh                                                # amount wiped off hands
      towasteb <- body1*removals$fwipeb                                                 # amount wiped off body
      hands2   <- hands1 - todrainh - towasteh                                          # amount left on hands after rinse/wipe
      body2    <- body1  - todrainb - towasteb                                          # amount left on body after rinse/wipe
      use.data[prod,3,y2=="fis",] <- use.data[prod,3,y2=="fis",] - todrains - towastes  # remove drain and waste from indoor surface 
      use.data[prod,3,y2=="fsk",] <- hands2 + body2                                     # amount left on skin
      todrain <- todrains + todrainh + todrainb                                         # total added to drain
      towaste <- towastes + towasteh + towasteb                                         # total added to waste
      use.data[prod,3,y2=="fdr",] <- use.data[prod,3,y2=="fdr",] + todrain              # add to existing drain amount
      use.data[prod,3,y2=="fws",] <- use.data[prod,3,y2=="fws",] + towaste              # add to existing waste amount
    }
    use.data[use.data<1E-9] <- 0                                                    # set all masses below 1E-9 milligrams to zero
    use.data[is.na(use.data)] <- 0                                                  # set undefined amounts to zero
    if (g$prog=="y") cat("\n  Evaluating use.data complete")
    if (g$save.r.objects=="y") {
      p   <- nrow(pucs)
      c   <- length(y2)
      num <- p*3*nc
      t   <- data.table(pucs$source.id[1],1,chem.list[1],matrix(0,nrow=num,ncol=c+1))
      setnames(t,c("puc","level","chem",y2,"tot"))
      row <- 0
      for (i in 1:nc) {
        for (k in 1:p) {
          for (j in 1:3) {
            row <- row+1
            t[row]$puc   <- pucs$source.id[k]
            t[row]$level <- j
            t[row]$chem  <- chem.list[i]
            t[row]$fia   <- use.data[k,j,y2=="fia",i]
            t[row]$fis   <- use.data[k,j,y2=="fis",i] 
            t[row]$fdr   <- use.data[k,j,y2=="fdr",i]
            t[row]$fws   <- use.data[k,j,y2=="fws",i] 
            t[row]$fsk   <- use.data[k,j,y2=="fsk",i]
            t[row]$foa   <- use.data[k,j,y2=="foa",i] 
            t[row]$fos   <- use.data[k,j,y2=="fos",i]
            t[row]$fap   <- use.data[k,j,y2=="fap",i] 
            t[row]$fsl   <- use.data[k,j,y2=="fsl",i]
            t[row]$fgi   <- use.data[k,j,y2=="fgi",i] 
            t[row]$fhn   <- use.data[k,j,y2=="fhn",i] 
            t[row]$tot   <- sum(use.data[k,j,1:c,i])
          }
        }
      }
      write.csv(t,paste0(g$out,"/Temp/use_data_",house.num,"_",g$run.name,".csv"),row.names=FALSE)
    }
    return(use.data)
  }
  
  
  
  # get.house.num extracts the house number from the file name
  
  get.house.num = function(houses,i) {
    h <- houses[i]
    a <- max(str_locate_all(h,"_")[[1]])+1
    b <- max(str_locate_all(h,".csv")[[1]])-4
    n <- as.numeric(substr(h,a,b))
    return(n)
  }
  
  
  
  # get.randoms returns random samples from a uniform distribution, using variable-specific seeds
  
  get.randoms = function(ran.vars,np,nc,all.chems,puc.types,flag) {
    # house variables (1 samples per house for each)
    hvar   <- ran.vars[[1]]
    nhvar  <- length(hvar)
    first  <- g$first.house
    last   <- g$last.house
    nsam   <- last-first+1
    hran   <- data.table(matrix(0,nrow=nsam,ncol=nhvar))
    setnames(hran,hvar)
    hseeds <- get.seeds(g$house.seed,nhvar)
    for (i in 1:nhvar) {
      b1 <- 2*i + flag-1
      b2 <- 2*i - flag
      set.seed(hseeds[b1:b2],"Marsaglia-Multicarry")
      hran[,i] <- runif(last)[first:last]
    }
    # puc variables (1 sample per house-puc combination, for each variable)
    pvar   <- unlist(unique(lapply(ran.vars[[2]],strip.n)))
    npvar  <- length(pvar)
    pran   <- data.table(matrix(0,nrow=nsam,ncol=npvar*np))
    setnames(pran,ran.vars[[2]])
    for (p in 1:np) {
      j      <- 953 * puc.types[source.id==puc.list[p]]$hem.id
      pseeds <- get.seeds(g$puc.seed+j,npvar)
      for (i in 1:npvar) {
        b1 <- 2*i + flag-1
        b2 <- 2*i - flag
        k  <- i + (p-1)*npvar
        y <- pseeds[b1:b2]
        set.seed(pseeds[b1:b2],"Marsaglia-Multicarry")
        pran[,k] <- runif(last)[first:last]
      }
    }
    # chem variables (1 sample per house-chem combination, for each variable)
    cvar   <- unlist(unique(lapply(ran.vars[[3]],strip.n)))
    ncvar  <- length(cvar)
    cran   <- data.table(matrix(0,nrow=nsam,ncol=ncvar*nc))
    setnames(cran,ran.vars[[3]])
    for (c in 1:nc) {
      j      <- 727 * all.chems[c]$chem.num
      cseeds <- get.seeds(g$chem.seed+j,ncvar)
      for (i in 1:ncvar) {
        b1 <- 2*i + flag-1
        b2 <- 2*i - flag
        k  <- i + (c-1)*ncvar
        set.seed(cseeds[b1:b2],"Marsaglia-Multicarry")
        cran[,k] <- runif(last)[first:last]
      }
    }
    house <- first:last
    return(data.table(house,hran,pran,cran))
  }
  
  
  
  # get.ranvars creates a list of variables requiring random draws
  
  get.ran.vars = function(nc,np) {
    if (g$prog=="y") cat("\n  Getting random vars...")
    ran.house      <- c("aer.out","height","lg.carb.f","lg.clean.air","lg.clean.sur","lg.depos","lg.load.air",
                        "lg.load.sur","lg.resus","sm.carb.f","sm.clean.air","sm.clean.sur","sm.depos",
                        "sm.load.air","sm.load.sur","sm.resus","temp","thick.bou","thick.sur","basal.vent",
                        "hand.wash","hand.mouth","thick.indir")
    prod.rv    <- c("prod.id","form.id","skin.frac")
    y <- ""
    for (i in 1:np) {
      y <- c(y,str_c(prod.rv,i))
    }  
    ran.prod <- y[-1]
    chem.rv        <- c("vapor","solub","kow","decay.air","decay.sur","diffus.air","conc.out")
    z <- ""
    for (j in 1:nc) {
      z <- c(z,str_c(chem.rv,j))
    }
    ran.chem <- z[-1]
    return(list(ran.house,ran.prod,ran.chem))
  }
  
  
  
  # get.seeds returns random number seeds for each variable
  
  get.seeds = function(init.seed,num) {
    n     <- 2*num
    base  <- as.integer64(2147483647)
    mult  <- as.integer64(397204094)
    seeds <- vector("integer",n)
    x     <- as.integer64(init.seed)
    for (i in 1:n) {
      x <- (x * mult) %% base
      seeds[i] <- as.numeric(x)
    }
    return(seeds)
  }
  
  
  
  # list.persons details all the persons in this household 
  
  list.persons = function(ph) {
    gens  <- vector("character",20)
    ages  <- vector("integer",20)
    prime <- rep(0,20)
    pset  <- 0
    pg   <- substr(ph$gender,1,1)
    pa   <- ph$age_years
    for (p in 1:20) {
      gens[p] <- substr(ph$genders,p,p)
      if(gens[p]==".") break
      ages[p] <- as.integer(substr(ph$ages,2*p-1,2*p))
      if(gens[p]==pg & ages[p]==pa & pset==0) {
        prime[p] <- 1
        pset     <- 1
      }  
    }
    np      <- p-1
    persons <- as.data.table(data.frame(1:np,gens[1:np],ages[1:np],prime[1:np]))
    setnames(persons,c("pnum","sex","age","primary"))
    persons$basal.vent[persons$primary==1] <- ph$basal.vent
    persons$skin.area[persons$primary==1]  <- ph$BSA_adj
    if(pset==0) cat("\n  Primary person not identified in this house") 
    return(persons)
  }
  
  
  
  # make.day.calendar creates a 364 day long standard calendar
  
  make.day.calendar = function() {
    c <- as.data.table(matrix(0,nrow=364,ncol=5))
    setnames(c,names(c),c("daynum","season","week","dayofweek","weekend"))
    c$daynum    <- 1:364
    c$season    <- 1+floor((c$daynum-1)/91)
    c$week      <- 1+floor((c$daynum-91*(c$season-1)-1)/7)
    c$dayofweek <- 1+((c$daynum-1) %% 7)
    c$weekend   <- 0
    c$weekend[c$dayofweek>5] <- 1
    return(c)  
  }
  
  
  
  # make.hour.calendar makes a calendar of 364 days of 24 hours each
  
  make.hour.calendar = function() {
    c <- as.data.table(matrix(0,nrow=8736,ncol=7))
    setnames(c,names(c),c("hournum","daynum","season","week","dayofweek","hourofday","weekend"))
    c$hournum   <- 1:8736
    c$daynum    <- 1+floor((c$hournum-1)/24)
    c$season    <- 1+floor((c$daynum-1)/91)
    c$week      <- 1+floor((c$daynum-91*(c$season-1)-1)/7)
    c$dayofweek <- 1+((c$daynum-1) %% 7)
    c$hourofday <- 1+((c$hournum-1) %% 24)
    c$weekend   <- 0
    c$weekend[c$dayofweek>5] <- 1
    return(c)  
  }
  
  
  
  # read.chem.fracs loads the .csv file of chemical fractions in each product
  
  read.chem.fracs = function(chem.list,puc.list,chem.frac.file) {
    if (is.null(chem.frac.file)) chem.frac.file <- g$chem.frac.file
    x <- fread(paste0("input/",chem.frac.file))
    setnames(x,tolower(names(x)))
    if(exists("id",x))  setnames(x,"id","source.id")
    if(exists("shedsid",x)) setnames(x,"shedsid","source.id")
    x$weight_fraction[is.na(x$weight_fraction)] <- 0
    if(length(puc.list)>0) x <- x[x$source.id %in% puc.list] 
    if(g$comp.method==2) x <- x[x$dtxsid %in% chem.list]
    if(g$comp.method==2) x <- x[x$weight_fraction>0]           # for comp.method=1, do this after making brand.list
    chem.fracs <- unique(x,by=c("product_id","formulation_id","dtxsid"))
    if (nrow(chem.fracs)==0) cat ("\n No product-chemical fractions left in run \n")
    return(chem.fracs)
  }
  
  
  
  
  # read.chem.props reads the .csv file of chemical properties 
  
  read.chem.props = function(chem.file,chem.list) {
    if (is.null(chem.file)) chem.file <- g$chem.file
    # Read chemical properties input file
    dt <- as.data.table(fread(paste0("input/",chem.file)))
    setnames(dt,names(dt),tolower(str_replace_all(names(dt),"_",".")))
    mode(dt$cas) <- "character"
    mode(dt$kp) <- "numeric"
    dt$cas <- trimzero(str_replace_all(dt$cas,"-","_"))
    dt$cas <- ifelse(substr(dt$cas,1,1)=="_",paste0("0",dt$cas),dt$cas)
    if(!exists("x",dt))               dt$x       <- 1:nrow(dt)
    if(is.na(dt[1]$x))                dt[1]$x    <- 1
    y <- lag(dt$x) + 1
    dt$x[is.na(dt$x)] <- y[is.na(dt$x)]
    if (length(chem.list)>0)  dt <- dt[dt$dtxsid %in% chem.list]
    dt <- unique(dt,by="dtxsid",fromLast=TRUE)
    if(nrow(dt)==0) cat("\n  No chemicals left in fugacity calcs \n")
    if(exists("name",dt))             setnames(dt,"name","chem.name")
    if(exists("vp.pa",dt))            setnames(dt,"vp.pa","vapor")
    if(exists("mw",dt))               setnames(dt,"mw","molwt")
    if(exists("water.sol.mg.l",dt))   dt$solub   <- dt$water.sol.mg.l/dt$molwt
    if(exists("log.kow",dt))          dt$kow     <- 10^dt$log.kow
    if(exists("half.sediment.hr",dt)) dt$decay.s <- 24*log(2)/dt$half.sediment.hr 
    if(exists("half.air.hr",dt))      dt$decay.a <- 24*log(2)/dt$half.air.hr 
    dt <- dt[order(dt$dtxsid,chem.list)]
    return(dt) 
  }
  
  
  
  # read.compart.fracs reads a .csv file of initial compartmental fractions by PUC
  
  read.compart.fracs = function(compart.file) {
    if (is.null(compart.file))  compart.file <- g$compart.file
    x <- fread(paste0("input/",compart.file))
    setnames(x,tolower(substr(names(x),1,3)))
    setnames(x,c("des","cod"),c("description","code"))
    x$code <- toupper(x$code)
    compart.fracs <- x
    return(compart.fracs)
  }
  
  
  
  # read.control file reads the settings for the current S2D run
  
  read.control.file = function(control.file) {
    if (is.null(control.file)) control.file <- "S2D_control_file.txt"
    lc <- str_length(control.file)
    if(!substr(control.file,lc-3,lc)=='.txt') control.file <- paste0(control.file,".txt")
    x <- fread(paste0("input/",control.file),header=FALSE,skip=0,col.names=c("key","setting"))
    x$key          <- tolower(x$key)
    chem.list      <- list(c(x$setting[x$key=="chem"]))
    puc.list       <- list(c(x$setting[x$key=="puc"]))
    fug.file       <- x$setting[x$key=="fug.file"] 
    chem.file      <- x$setting[x$key=="chem.file"]
    chem.frac.file <- x$setting[x$key=="chem.frac.file"]
    puc.type.file  <- x$setting[x$key=="puc.type.file"]
    compart.file   <- x$setting[x$key=="compart.file"]
    puc.met.file   <- x$setting[x$key=="puc.met.file"]
    skin.area.file <- x$setting[x$key=="skin.area.file"]
    removal.file   <- x$setting[x$key=="removal.file"]
    vent.file      <- x$setting[x$key=="vent.file"]
    diary.prefix   <- x$setting[x$key=="diary.prefix"]
    run.name       <- x$setting[x$key=="run.name"]
    house.seed     <- x$setting[x$key=="house.seed"]
    puc.seed       <- x$setting[x$key=="puc.seed"]
    chem.seed      <- x$setting[x$key=="chem.seed"]
    init.seed      <- x$setting[x$key=="init.seed"]
    if (length(init.seed)>0) { 
      house.seed   <- init.seed
      puc.seed     <- init.seed
      chem.seed    <- init.seed
    } else init.seed <- 0
    mode(house.seed)  <- "integer"
    mode(puc.seed)    <- "integer"
    mode(chem.seed)   <- "integer"
    first.house       <- x$setting[x$key=="first.house"]
    if (length(first.house)==0) first.house <- 1
    mode(first.house) <- "integer"
    last.house        <- x$setting[x$key=="last.house"]
    if (length(last.house)==0) last.house <- NA
    mode(last.house)  <- "integer"
    n.houses          <- x$setting[x$key=="n.houses"]
    if (length(n.houses)==0)   n.houses   <- NA
    mode(n.houses)    <- "integer"
    if (is.na(last.house) & !is.na(n.houses)) last.house <- first.house+n.houses-1
    if (is.na(n.houses) & !is.na(last.house)) n.houses   <- last.house-first.house+1
    comp.method       <- x$setting[x$key=="comp.method"]
    if (length(comp.method)==0) comp.method <- 1
    mode(comp.method) <- "integer" 
    puc.offset        <- x$setting[x$key=="puc.offset"]
    if (length(puc.offset)==0) puc.offset <- 0
    mode(puc.offset)  <- "integer"
    prog              <- substr(x$setting[x$key=="show.progress"],1,1)
    if (length(prog)==0)  prog <- "n"                     # default is no progress messages
    parallel          <- substr(x$setting[x$key=="parallel"],1,1)
    if (length(parallel)==0) parallel <- "n"              # default is not parallel
    save.r.objects    <- substr(x$setting[x$key=="save.r.objects"],1,1)
    if (length(save.r.objects)==0) save.r.objects <- "n"  # default is not to save
    out <- paste0("output/S2D/",run.name)
    g <- as.data.table(list(chem.list,puc.list,fug.file,chem.file,chem.frac.file,puc.type.file,compart.file,puc.met.file,
                            skin.area.file,removal.file,vent.file,diary.prefix,run.name,n.houses,init.seed,prog,parallel,
                            puc.offset,first.house,last.house,comp.method,save.r.objects,house.seed,puc.seed,chem.seed,out))
    setnames(g,c("chem.list","puc.list","fug.file","chem.file","chem.frac.file","puc.type.file","compart.file",
                 "puc.met.file","skin.area.file","removal.file","vent.file","diary.prefix","run.name","n.houses",
                 "init.seed","prog","parallel","puc.offset","first.house","last.house","comp.method","save.r.objects",
                 "house.seed","puc.seed","chem.seed","out"))
    return(g)
  }
  
  
  
  # read.diary loads the product-use diary for one household
  
  read.diary = function(diary.prefix,run.name,house.num,persons) {
    if (is.null(diary.prefix)) diary.prefix <- g$diary.prefix
    if (is.null(run.name))         run.name <- g$run.name
    diary <- fread(paste0("output/ABM/", diary.prefix, house.num, ".csv"),stringsAsFactors = FALSE, na.strings = c("","NA"))
    setnames(diary,tolower(names(diary)))
    if(exists("person.index",diary))                        setnames(diary,"person.index","p")
    if(exists("day.of.the.year",diary))                     setnames(diary,"day.of.the.year","daynum")
    if(exists("sheds.id.refined",diary))                    setnames(diary,"sheds.id.refined","source.id")
    if(exists("start.time.hr.using.military.time",diary))   setnames(diary,"start.time.hr.using.military.time","start")
    if(exists("use.ht",diary))                              setnames(diary,"use.ht","hand.dur")
    if(exists("use.mass",diary))                            setnames(diary,"use.mass","mass")
    if(exists("household.index",diary))                     setnames(diary,"household.index","house.num")
    if(exists("person.gender",diary))                       setnames(diary,"person.gender","sex")
    if(exists("person.age",diary))                          setnames(diary,"person.age","age")
    if(exists("primary.person",diary))                      setnames(diary,"primary.person","primary")
    if(exists("indoor.outdoor",diary))                      setnames(diary,"indoor.outdoor","indoor_outdoor")
    if(exists("pucid.productype",diary))                    setnames(diary,"pucid.productype","product_type")
    primenum <- persons$pnum[persons$primary==1]
    diary$person[diary$primary==1] <- primenum
    diary$person[diary$primary==0 & diary$p<=primenum] <- diary$p[diary$primary==0 & diary$p<=primenum]-1
    diary$person[diary$primary==0 & diary$p >primenum] <- diary$p[diary$primary==0 & diary$p >primenum]           
    diary$mass[is.na(diary$mass)] <- 0
    diary$source.id[is.na(diary$source.id)] <- "none"
    diary$mass[diary$source.id=="none"]     <- 0
    diary$hand.dur[diary$source.id=="none"] <- 0
    diary$product_type[diary$source.id=="none"] <- "none"
    diary$indoor_outdoor[diary$source.id=="none"] <- "I"
    diary$hour           <- 0
    diary$hour <- 1 + floor(diary$start)
    diary$hournum  <- 24*(diary$daynum-1) + diary$hour
    setorder(diary,person,hournum,start)
    diary$row <- 1:nrow(diary)
    return(diary)
  }
  
  
  
  # read.fug.inputs reads the .csv file of chemical-independent fugacity variables
  
  read.fug.inputs = function(fug.file) {
    if (is.null(fug.file)) fug.file <- g$fug.file
    # Read fugacity input file
    df <- read.csv(paste0("input/",fug.file),as.is=TRUE) 
    names(df) <- tolower(names(df))
    if (nrow(df)==0) stop ("No data on fugacity input file \n")
    options(warn=-1)
    mode(df$par1) <- "numeric"
    mode(df$par2) <- "numeric" 
    mode(df$par3) <- "numeric"
    mode(df$par4) <- "numeric"
    mode(df$lower.trun) <- "numeric"
    mode(df$upper.trun) <- "numeric"
    options(warn=0)
    df <- transform(df, form = tolower(gsub(" ", "", form)))
    dt <- as.data.table(df)
    dt$varname <- tolower(dt$varname)
    return(dt)
  }
  
  
  
  # read.pophouse reads the output from the HEM RPGen module 
  
  
  read.pophouse = function(run.name) {
    if (g$prog=="y") cat("\n  Reading pophouse...")
    if (is.null(run.name)) run.name <- g$run.name
    pophouse <- fread(paste0("output/RPGen/pophouse.csv"))
    if (!exists("basa.vent",pophouse)) pophouse$basal.vent <- 5+0.5*pmin(14,floor(pophouse$age_years/2))
    pophouse$unitsf[pophouse$unitsf<0] <- 1500
    pophouse$house <- 1:nrow(pophouse)
    return(pophouse)
  }
  
  
  
  # read.puc.met reads the input file containing the mean MET value for each PUC
  
  read.puc.met = function(puc.met.file,puc.list) {
    if (is.null(puc.met.file)) puc.met.file <- g$puc.met.file
    x <- fread(paste0("input/",puc.met.file))
    setnames(x,tolower(names(x)))
    if (exists("new_puc",x)) setnames(x,"new_puc","source.id")
    if(length(puc.list>0)) {
      x <- x[x$source.id %in% puc.list] 
      for (i in 1:length(puc.list)) {
        if (!puc.list[i] %in% x$source.id) cat("  PUC ",puc.list[i]," not on PUC MET file")
      }
    }
    setnames(x,c("new_general_category","new_product type","new_refined product type","chad activity","mean mets from apex"),
             c("category","product_type","refined","chad_act","met"))
    return(x)
  }
  
  
  
  # read.puc.types reads the .csv file listing the product type for each PUC
  
  read.puc.types = function(puc.type.file,puc.list) {
    if (is.null(puc.type.file)) puc.type.file <- g$puc.type.file
    x <- fread(paste0("input/",puc.type.file))
    setnames(x,tolower(names(x)))
    if (exists("new_source.id",x)) setnames(x,"new_source.id","source.id")
    if(length(puc.list>0)) {
      x <- x[x$source.id %in% puc.list] 
      for (i in 1:length(puc.list)) {
        if (!puc.list[i] %in% x$source.id) cat("  PUC ",puc.list[i]," not on PUC types file")
      }
    }
    x$code    <- toupper(x$code)
    puc.types <- x[x$code!="XXX"]
    setorder(puc.types,source.id)
    if (g$save.r.objects=="y") {
      write.csv(puc.types,paste0(g$out,"/Temp/puc_types_",run.name,".csv"))
    }
    return(puc.types)
  }
  
  
  
  #read.removals reads the file with rinse-off and wipe-off fractions by PUC
  
  read.removals = function(removal.file) {
    if (is.null(removal.file)) removal.file <- g$removal.file
    x <- fread(paste0("input/",removal.file))
    setnames(x,tolower(names(x)))
    if (exists("new_puc",x)) setnames(x,"new_puc","source.id")
    return(x)
  }
  
  
  
  # read.skin.areas reads the file with affected skin area fractions by PUC
  
  read.skin.areas = function(skin.area.file) {
    if (is.null(skin.area.file)) skin.area.file <- g$skin.area.file  
    x <- fread(paste0("input/",skin.area.file))
    setnames(x,tolower(names(x)))
    if (exists("new_puc",x)) setnames(x,"new_puc","source.id")
    skin.areas <- unique(x,by="source.id")
    return(skin.areas)
  }
  
  
  # read.vent.file reads the .csv file with basal ventilation rates in (m3/day)
  
  read.vent.file = function(vent.file){
    if (is.null(vent.file)) vent.file <- g$vent.file
    x <- fread(paste0("input/",vent.file))
    setnames(x,tolower(names(x)))
    return(x)
  }
  
  
  
  # shorten brand list if comp.method==2
  
  reeval.brand.list = function(brand.list,chem.list,chem.fracs) {
    for (i in 1:nrow(brand.list)) {
      x <- brand.list[i]
      x$brands[[1]] <- list(unique(chem.fracs$product_id[chem.fracs$source.id==brand.list[i]$puc]))
      if(i==1) {y <- x} else {y <- rbind(y,x)}
    }
    return(y)
  }
  
  
  
  # string.values converts a string containing numbers to a list of numbers
  
  string.values = function(x){
    w <- lapply(str_split(x," "),as.numeric)[[1]]
    v <- w[!is.na(w)]
    return(v)
  }
  
  
  
  # strip.n removes any numerical suffix from a variable name
  
  strip.n = function(string) {
    new <- string
    y <- substr(new,str_length(new),str_length(new))
    while (y %in% c("0","1","2","3","4","5","6","7","8","9")) {
      new <- substr(new,1,str_length(new)-1)
      y <- substr(new,str_length(new),str_length(new))
    }
    return(new)
  }
  
  
  
  # trimzero removes leading zeroes from CAS numbers
  
  trimzero = function(x) {
    y <- x
    for (j in 1:length(x)) {
      i <- 1
      while(substr(x[j],i,i)=="0" & substr(x[j],i+1,i+1)!="_") {
        substr(y[j],i,i)<- " "
        i <- i+1
      }  
    }
    return(str_trim(y))
  }
  
  
  
  write.all.annual = function() {
    ann.list <- list.files(paste0(g$out,"/Annual"),pattern="^[House_]")
    n <- length(ann.list)
    if (n>0) all.annual <- as.data.table(read.csv(paste0(g$out,"/Annual/",ann.list[1])))
    if (n>1) {
      for (i in 2:length(ann.list)) {
        x    <- as.data.table(read.csv(paste0(g$out,"/Annual/",ann.list[i])))
        all.annual <- rbind(all.annual,x)
      }
    }
    if(exists("all.annual")) {
      setorder(all.annual,household)
      write.csv(all.annual,paste0(g$out,"/Annual/All_houses_annual.csv"),row.names=FALSE)
    }  
    return(cat("\nS2D finished\n\n"))
  }   
  
  
  
  write.all.houses = function() {
    houses <- list.files(paste0(g$out,"/Chem"),pattern="^[House_]")
    if (length(houses)>0) {
      all.houses    <- fread(paste0(g$out,"/Chem/",houses[1]))
      all.houses$V1 <- get.house.num(houses,1)
    }
    if (length(houses)>1) {
      for (i in 2:length(houses)) {
        x    <- fread(paste0(g$out,"/Chem/",houses[i]))
        x$V1 <- get.house.num(houses,i)
        all.houses <- rbind(all.houses,x)
      }
    }
    if(exists("all.houses")) {
      setnames(all.houses,1,"house.num")
      setorder(all.houses,house.num)
      write.csv(all.houses,paste0(g$out,"/Chem/All_houses.csv"),row.names=FALSE)
    }  
  }
  
  
  
  
  
  
  
  
  
  S2D_Runner = function(i) {
    library(data.table)
    library(stringr)
    library(plyr)
    library(dplyr)
    library(dtplyr)
    library(ggplot2)
    library(bit64)
    library(foreach)
    library(doParallel)
    library(reshape2)
    
    house.num          <- i
    if (g$prog=="y") cat("\nHouse=", house.num)
    pd                 <- person.data[person.data$house==house.num]
    hp                 <- house.props[house.props$house==house.num]
    q                  <- q.randoms[q.randoms$house==house.num]
    persons            <- list.persons(pd)
    prime              <- eval.prime(persons,ventilation,q)
    diary              <- read.diary(g$diary.prefix,g$run.name,house.num,persons)
    diary.pucs         <- distinct(diary,source.id,.keep_all=TRUE)
    diary.pucs         <- select(diary.pucs,one_of(c("source.id","product_type","indoor_outdoor","hand.dur","mass")))
    setorder(diary.pucs,source.id)
    pucs               <- left_join(diary.pucs[diary.pucs$source.id %in% puc.list],puc.codes,by="source.id")
    mode(pucs$mass)    <- "numeric"
    direct             <- 0
    indirect           <- 0
    env.impact         <- 0
    lcia               <- NULL
    chem.props         <- eval.chem.props(fug.cvars,chem.list,ran.vars,q)
    flows              <- eval.flows(hp,chem.props)
    setorder(pucs,source.id)
    pucs$met           <- puc.codes$met[puc.types$source.id %in% pucs$source.id]
    pucs$spray         <- substr(pucs$code,1,1)=="S"
    pucs.areas         <- eval.pucs.areas(pucs,skin.areas,q)
    puc.wipe.rinse     <- eval.puc.wipe.rinse(pucs,removal.fracs)
    prod.chem          <- eval.prod.chem(pucs,puc.list,brand.list,chem.list,chem.fracs,release.fracs,q,house.num)
    dermal.rates       <- eval.dermal.rates(fug.cvars,hp,prime,pucs.areas,prod.chem,nc,house.num)
    chem.release       <- NULL
    indoor.gaps        <- NULL
    indoor.hrs         <- NULL
    use.data           <- NULL
    d                  <- diary[diary$person==prime$pnum]
    chem.totals        <- 0
    if (nrow(pucs)>0) {
      use.chem           <- eval.use.chem(pucs,prod.chem,nc)
      if (sum(use.chem)>0) {
        compart.list     <- select_vars(names(prod.chem),starts_with("f"),exclude="formulation_id")
        use.data         <- eval.use.data(pucs,prod.chem,fug.cvars,dermal.rates,puc.wipe.rinse,nc,house.num)
        chem.release     <- eval.chem.release(pucs,use.data,use.chem,compart.list,diary,nc,calendar.hours)
        if (g$save.r.objects=="y" & any(chem.release[8:ncol(chem.release)]>0)) {
          write.csv(chem.release,paste0(g$out,"/Temp/chem.release_",house.num,"_",run.name,".csv"))
        } 
        chem.totals      <- eval.chem.totals(chem.release,chem.list,house.num)
        release.times    <- eval.release.times(chem.release)
        release.hrs      <- unlist(release.times[1])
        release.gaps     <- unlist(release.times[2])
        indoor.hrs       <- unlist(release.times[3])
        indoor.gaps      <- unlist(release.times[4])
        direct           <- eval.direct(d,use.data,use.chem,pucs,dermal.rates,compart.list,prime,puc.wipe.rinse,chem.list,fug.cvars,chem.totals)
        if (g$save.r.objects=="y" & any(direct[2:ncol(direct)]>0)) write.csv(direct,paste0(g$out,"/Temp/direct_",house.num,"_",run.name,".csv"))
      }
    }  
    fug.day        <- as.data.table(matrix(0,nrow=364,ncol=4*nc+2))
    setnames(fug.day,c("daynum","hour",str_c(c(rep("air",nc),rep("sur",nc),rep("win",nc),rep("was",nc)),rep(1:nc,4))))
    fug.day$daynum <- 1:364
    fugs           <- eval.fug.concs.an(chem.release,flows,indoor.hrs,indoor.gaps)
    fug.day        <- select(fugs[1:364],c(1,3:(4*nc+2)))
    if (g$save.r.objects=="y" & any(fug.day[,2:ncol(fug.day)]>0)) write.csv(fug.day,paste0(g$out,"/Temp/fug.day_",house.num,"_",run.name,".csv"))
    fug.hour       <- select(fugs[365:9100],c(1:(2*nc+2)))
    if (g$save.r.objects=="y" & any(fug.hour[,2:ncol(fug.hour)]>0)) write.csv(fug.hour,paste0(g$out,"/Temp/fug.hour_",house.num,"_",run.name,".csv"))
    indirect       <- eval.indirect(d,fug.hour,dermal.rates,hp,nc,prime,fug.cvars,chem.list)
    if (g$save.r.objects=="y" & any(indirect[2:ncol(indirect)]>0)) write.csv(indirect,paste0(g$out,"/Temp/indirect_",house.num,"_",run.name,".csv"))
    env.impact       <- eval.env.impact(diary,use.data,pucs,fug.day,compart.list,chem.totals,nc)
    if (g$save.r.objects=="y" & ncol(env.impact)>1 & any(env.impact[2:ncol(env.impact)]>0)) {
      write.csv(env.impact,paste0(g$out,"/Temp/env.impact_",house.num,"_",run.name,".csv"))
    }
    annual      <- eval.summary(direct,indirect,env.impact,g$run.name,house.num,nc,chem.list,chem.totals)
    if (g$save.r.objects=="y" & any(as.data.frame(annual)[3:ncol(annual)]>0)) write.csv(annual,paste0(g$out,"/Temp/annual_",house.num,"_",run.name,".csv"))
    if (length(puc.list)==1 & nrow(pucs)==1) {
      lcia <- eval.lcia(puc.list,chem.list,prod.chem,diary,annual,house.num,nc)
      if (any(lcia$chem.mass>0)) write.csv(lcia,file=paste0(g$out,"/LCIA/LCIA_",house.num,".csv"),row.names=FALSE)
    }  
    return(lcia)
  }    
  
  
  ###################################################################
  ############## End of function definitions ########################
  ###################################################################
  
  
  
  
  
  
  
  
  
  
  ###################################################################
  ############## Start of S2D executable statements #################
  ###################################################################
  
  
  
  #S2D preprocessing code before looping on households
  g             <- read.control.file(control.file)
  run.name      <- g$run.name
  if(!dir.exists("output/S2D"))               dir.create("output/S2D")
  if(!dir.exists(g$out))                      dir.create(g$out)
  if(!dir.exists(paste0(g$out,"/Daily")))     dir.create(paste0(g$out,"/Daily"))
  if(!dir.exists(paste0(g$out,"/Annual")))    dir.create(paste0(g$out,"/Annual"))
  if(!dir.exists(paste0(g$out,"/LCIA")))      dir.create(paste0(g$out,"/LCIA"))
  if(!dir.exists(paste0(g$out,"/Prod_chem"))) dir.create(paste0(g$out,"/Prod_chem"))
  if(!dir.exists(paste0(g$out,"/Chem")))      dir.create(paste0(g$out,"/Chem"))
  if(!dir.exists(paste0(g$out,"/Temp")))      dir.create(paste0(g$out,"/Temp"))
  if (g$prog=="y") cat("\nRun name = ",run.name)  
  if(exists("number.of.houses")) { if(!is.null(number.of.houses)) g$n.houses <- number.of.houses }
  house.list     <- g$first.house:g$last.house
  chem.list      <- unlist(g$chem.list)
  puc.list       <- unlist(g$puc.list)
  if (length(chem.list)==0) chem.list <- NULL
  if (length(puc.list)==0)   puc.list <- NULL
  if (g$prog=="y") cat("\nReading data...")  
  fug.hvars      <- read.fug.inputs(g$fug.file) 
  fug.cvars      <- read.chem.props(g$chem.file,chem.list)
  if(length(chem.list)==0 | is.null(chem.list)) chem.list <- unique(fug.cvars$dtxsid)
  
  puc.types      <- read.puc.types(g$puc.type.file,puc.list)
  puc.met        <- read.puc.met(g$puc.met.file,puc.list)
  puc.codes      <- join(select(puc.types,source.id,code),select(puc.met,source.id,met),by="source.id")
  puc.list       <- as.vector(puc.types$source.id)
  chem.fracs     <- read.chem.fracs(chem.list,puc.list,g$chem.frac.file)
  brand.list     <- eval.brand.list(chem.list,puc.list,chem.fracs)
  puc.list       <- brand.list$puc
  chem.fracs     <- chem.fracs[chem.fracs$source.id %in% puc.list]
  chem.list      <- eval.chem.list(chem.list,puc.list,chem.fracs)
  if (length(chem.list)==0) {
    stop("\n No chemicals left to model\n")
  }
  chem.fracs     <- chem.fracs[chem.fracs$dtxsid %in% chem.list]
  if(g$comp.method==2) brand.list <- reeval.brand.list(brand.list,chem.list,chem.fracs)
  all.chems      <- data.table(1:length(chem.list),chem.list)
  setnames(all.chems,c("chem.num","dtxsid"))
  all.chems      <- inner_join(select(fug.cvars,dtxsid,x),all.chems,by="dtxsid")
  setorder(all.chems,chem.num)
  fug.cvars      <- fug.cvars[fug.cvars$dtxsid %in% chem.list]
  fug.cvars      <- fug.cvars[order(fug.cvars$dtxsid,chem.list)]
  
  ventilation    <- read.vent.file(g$vent.file)
  skin.areas     <- read.skin.areas(g$skin.area.file)
  removal.fracs  <- read.removals(g$removal.file)
  compart.fracs  <- read.compart.fracs(g$compart.file)
  release.fracs  <- eval.release.fracs(puc.types,compart.fracs)
  pophouse       <- read.pophouse(g$run.name)
  if(is.na(g$n.houses))   g$n.houses   <- nrow(pophouse)
  if(is.na(g$last.house)) g$last.house <- g$n.houses
  person.data    <- pophouse[g$first.house:g$last.house]
  if (g$prog=="y") cat("\nDone with reading input data...")
  
  nc             <- length(chem.list)
  np             <- length(puc.list)
  ran.vars       <- get.ran.vars(nc,np)
  q.randoms      <- get.randoms(ran.vars,np,nc,all.chems,puc.types,flag=0)
  if(g$save.r.objects=="y") write.csv(q.randoms,paste0(g$out,"/Temp/q.randoms_",run.name,".csv"))
  house.props    <- eval.house.props(fug.hvars,person.data,q.randoms)
  calendar.hours <- make.hour.calendar()
  calendar.days  <- make.day.calendar()
  lcia.count     <- 0
  cols           <- 23*nc+1
  all.houses     <- data.table(matrix(0,nrow=g$n.houses,ncol=cols))
  all.lcia <- NULL
  if (g$prog=="y") cat("\nStarting ",g$n.houses," houses...")
  
  
  # Call S2D runner
  if (g$parallel=="y") {
    # register number of cores (total #CPU-2)
    cl <- makeCluster(detectCores() - 2)
    registerDoParallel(cl, cores = detectCores() - 2)
    writeLines(c(""), "log.txt")
    out2 <- foreach(i = g$first.house:g$last.house) %dopar% {
      sink(paste0("output/S2D/log",i,".txt"), append=FALSE)  
      res <- tryCatch({ 
        cat(paste("Starting iteration",i,"\n"))  
        S2D_Runner(i) 
      }, error=function(e) NULL)
      sink()
    }
    on.exit(stopCluster(cl))
    all.lcia <- rbindlist(out2)
  }
  
  
  if (g$parallel=="n") {
    for (i in g$first.house:g$last.house) {
      lcia <- S2D_Runner(i)
      if (exists("lcia")) {
        if (!exists("all.lcia")) all.lcia <- lcia else all.lcia <- rbind(all.lcia,lcia)
      }  
    } 
  }  
  ########S2D post processing to get lcia.avg after all households are run
  if (length(puc.list)==1 & exists("all.lcia") & length(all.lcia)>0) {
    lcia.avg <- eval.lcia.avg(all.lcia,puc.list,chem.list,nc)
    write.csv(lcia.avg,file=paste0("output/S2D/",g$run.name,"/LCIA/LCIA_avg.csv"),row.names=FALSE)
  }
  write.all.houses()
  write.all.annual()
}

###################################################################
############## End of S2D code ####################################
###################################################################

