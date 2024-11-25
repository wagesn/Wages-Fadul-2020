
###install required R packages
library(Iso)
library(nnet)

###Load the function 'WFdesign' 
WFdesign<-function(truth,fprob,a0,b0,q.skel,aplus,tul,ell,Mmax,Nmax,start,put,puf){

### run a trial 	
    ndose = length(truth);   #number of combos
    y=n=z=M=dose.select=ptox.hat=peff.hat=p.infeasible=rep(0,ndose); 
    jstar = start;   #current dose level	 
    stopf=stops=untx=0; #indicate if trial stops early
    i=1	
    A=q.skel*aplus
    B=aplus-A

while(i <= Mmax)
    {
	L=runif(1)
	fdoses=which(fprob>L)
	Y=ifelse(length(fdoses)==0,0,max(fdoses))
	if(Y==0){untx=untx+1}
	z[fdoses]=z[fdoses]+1
	M=M+1
	
	for(j in 1:ndose){
		p.infeasible[j] = pbeta(ell, z[j] + A[j], M[j] - z[j] + B[j]);
        }

	if(p.infeasible[1]> puf){
			stopf=1
			break
		}


if(Y>0){
	curr=min(jstar,Y)
	y[curr] = y[curr] + rbinom(1,1,truth[curr]);
	n[curr] = n[curr] + 1;
	tried=which(n>0)

	if(1 - pbeta(tul, y[1] + a0, n[1] - y[1] + b0) > put){
			stops=1
			break
		}

		if(length(tried)<ndose & all(y==0)){
			jstar<-jstar+1
		} else {
			u=(y[tried]+a0)/(n[tried]+a0+b0)
			pipost=pava(u,w=n[tried])
			lossvec=ifelse(pipost > tul, (1-0.5)*(pipost-tul), 0.5*(tul-pipost))  
			T=lossvec==min(lossvec)
			poss=which(T)
			if(sum(T)==1){
				sugglev=poss
			} else {
				if(all(pipost[poss]>tul)){
					sugglev=min(poss)
				} else {
					sugglev=max(poss)
				}
			}
			if(pipost[sugglev]<tul & length(tried)<ndose){
				jstar=ifelse(n[sugglev+1]==0,sugglev+1,sugglev)					
			} else {
				jstar=sugglev
			}
		}
	}

	fset=which(p.infeasible<puf)
	jstar=min(max(fset),jstar)

	if(sum(n)>=Nmax){
		stops=0
		break
		}
	i=i+1
}
	if(stops==0 & stopf==0){
		fmtd=min(max(fset),jstar)
		dose.select[fmtd]=dose.select[fmtd]+1;
		}
	return(list(dose.select=dose.select,tox.data=y,pt.allocation=n,stopf=stopf,stops=stops,unt=untx))
}
##########'WFdesign' end here

###Load the function 'WFdesign.sim' 
WFdesign.sim<-function(ntrial,truth,fprob,a0,b0,q.skel,aplus,tul,ell,Mmax,Nmax,start,put,puf){
	ndose=length(truth)
	nuntx=rep(0,ntrial)
	dose.select<-y<-z<-n<-naf<-matrix(nrow=ntrial,ncol=ndose)
	nstopf=nstops=0
	
	for(i in 1:ntrial){
		result<-WFdesign(truth,fprob,a0,b0,q.skel,aplus,tul,ell,Mmax,Nmax,start,put,puf)
		dose.select[i,]=result$dose.select
		y[i,]=result$tox.data
		n[i,]=result$pt.allocation
		#z[i,]=result$feas.data
		#naf[i,]=result$feas.eval
		nuntx[i]=result$unt
		nstopf=nstopf+result$stopf
		nstops=nstops+result$stops
	}
	cat("True tox probability:\n");
	cat(round(truth,3), sep="\t",  "\n");
	cat("True feas probability:\n");
	cat(round(fprob,3), sep="\t",  "\n");
	cat("FMTD selection percentage:\n");
	cat(formatC(colMeans(dose.select)*100, digits=1, format="f"), sep="\t",  "\n");
	cat("Average nmber of DLTs:\n");
	cat(formatC(colMeans(y), digits=1, format="f"), sep="\t",   "\n");
	cat("Average number of patients infused:\n");
	cat(formatC(colMeans(n), digits=1, format="f"), sep="\t",   "\n");
	#cat("number feasible:   		    			", formatC(colMeans(z), digits=1, format="f"), sep="\t",   "\n");
	#cat("number evaluated feasibility:      			", formatC(colMeans(M), digits=1, format="f"), sep="\t",   "\n");
	cat("Average number of patients not treated:\n");
	cat(formatC(mean(nuntx), digits=1, format="f"), sep="\t",   "\n");
	cat("percentage of stop (safety):\n");
	cat(nstops/ntrial*100, "\n");
	cat("percentage of stop (feasibility):\n");
	cat(nstopf/ntrial*100, "\n");
} 
##########'WFdesign.sim' end here


##Comparison to Thall et al. (2001)
t1<-c(0.10,0.30,0.50,0.70,0.80)
f1<-c(0.99,0.95,0.90,0.75,0.50)

t2<-c(0.10,0.30,0.50,0.70,0.80)
f2<-c(0.90,0.75,0.50,0.25,0.05)

t3<-c(0.05,0.10,0.30,0.50,0.60)
f3<-c(0.25,0.10,0.05,0.02,0.01)

t4<-c(0.01,0.05,0.07,0.10,0.30)
f4<-c(0.99,0.95,0.90,0.75,0.50)

t5<-c(0.50,0.60,0.70,0.75,0.80)
f5<-c(0.90,0.75,0.50,0.25,0.05)

tul=0.3            ##target toxicity rate 
ell=0.5		 ##minimum infusibility threshold
Nmax=24            ##maximum number evaluated for toxicity
Mmax=48		 ##maximum total sample size
start=1            ##starting combination
ntrial=10000       ##number of simulated trials 
put=0.7		 ##upper probability cutoff for safety
puf=0.7            ##upper probability cutoff for feasibility
q.skel<-c(0.90,0.85,0.80,0.75,0.70) ##prior means for infusibility probabilities
aplus=rep(0,length(q.skel))	      ##prior sample size at each dose level

###calculate prior for toxicity
x<-2*tul
mu<-tul
u<-0.95
f<-function(b){
	pbeta(x,mu*b/(1-mu),b)-u
}
b0<-uniroot(f,c(0.0001,100))$root
a0<-mu*b0/(1-mu)
round(c(a0,b0),2)

###calculate prior for unfusibility
for(i in 1:length(q.skel)){
	x=q.skel[i]/2
	mu=q.skel[i]
	f<-function(b){
		1-pbeta(x,mu*b/(1-mu),b)-0.95
	}
	bF=uniroot(f,c(0.0001,100))$root
	aF=mu*bF/(1-mu)
aplus[i]=round(aF+bF,2)
}

##simulate a single trial
#WFdesign(truth,fprob,a0,b0,q.skel,aplus,tul,ell,Mmax,Nmax,start,put,puf)

truth<-t2 ##true toxicity probability scenario
fprob<-f2 ##true feasibility probability scenario
WFdesign.sim(ntrial,truth,fprob,a0,b0,q.skel,aplus,tul,ell,Mmax,Nmax,start,put,puf)


##Application to the glioblastoma phase I trial 
d<-4
###True toxicity probability scenarios
r1<-c(0.01,0.04,0.07,0.12)
r2<-c(0.04,0.07,0.12,0.25)
r3<-c(0.07,0.12,0.25,0.36)
r4<-c(0.12,0.25,0.36,0.50)
r5<-c(0.25,0.36,0.50,0.65)
r6<-c(0.50,0.65,0.79,0.85)

###True un-feasibility probability scenarios
nf1<-rep(1,d)
nf1<-c(0.99,0.97,0.95,0.90)
nf2<-c(0.95,0.90,0.80,0.70)
nf3<-c(0.90,0.85,0.68,0.60)
nf4<-c(0.85,0.70,0.65,0.55)
nf5<-c(0.65,0.55,0.45,0.35)

tul=0.25            ##target toxicity rate 
ell=0.8		  ##minimum infusibility threshold
Nmax=24             ##max number of pts evaluated for toxicity 	
Mmax=30		  ##max total sample size
start=1             ##starting combination
ntrial=1000        ##number of simulated trials 
put=0.9		  ##upper probability cutoff for safety
puf=0.7             ##upper probability cutoff for feasibility 
q.skel<-c(0.90,0.85,0.80,0.75)	##prior means for infusibility probabilities
aplus=rep(0,length(q.skel))	      ##prior sample size at each dose level

###calculate prior for unfusibility
for(i in 1:length(q.skel)){
	x=q.skel[i]/2
	mu=q.skel[i]
	f<-function(b){
		1-pbeta(x,mu*b/(1-mu),b)-0.95
	}
	bF=uniroot(f,c(0.0001,100))$root
	aF=mu*bF/(1-mu)
aplus[i]=round(aF+bF,2)
}

###calculate prior for toxicity
x<-2*tul
mu<-tul
u<-0.95
f<-function(b){
	pbeta(x,mu*b/(1-mu),b)-u
}
b0<-uniroot(f,c(0.0001,100))$root
a0<-mu*b0/(1-mu)
round(c(a0,b0),2)

##simulate a single trial
#WFdesign(truth,fprob,a0,b0,q.skel,aplus,tul,ell,Mmax,Nmax,start,put,puf)

truth<-r1  ##true toxicity probability scenario
fprob<-nf1 ##true feasibility probability scenario
WFdesign.sim(ntrial,truth,fprob,a0,b0,q.skel,aplus,tul,ell,Mmax,Nmax,start,put,puf)


