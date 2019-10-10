// SSM with priors for Pacific Ocean Perch (POP)
// cf. Appendix F (p. 76) of Edwards, Haigh and Starr (2014) DFO CSAS Res Doc 2013/093
// v0.4
#include <TMB.hpp>


/*
### Indices
  - a=1,...,A: age class, A=30
  - t=1,...,TC: time, TC=73, corresponds to years 1940,...,2012.
    t=0 only for B0 (separate from Bt) initial unfished equilibrium conditions
  - g=1,2,3,G=4: fleet/data source:
		  1=commercial trawl data
      2=West Coast Vancouver Island synoptic survey series
      3=National Marine Fisheries Service Triennial survey series
      4=GB Reed historical survey series
  - s: sex, 1=females, 2=males
### Random effects summary
  - R_t is the only random effect
  - N_{a,t,s} is a deterministic function of C_t, w_{a,s}, s_{a,g,s}, R_0, M_s and R_t
  - B_t is a deterministic function of w_{a,s}, m_a and predicted N_{a,t,s}
### Fixed param summary
  - s_{a,g,s} is a deterministic function of mu_g, delta_g and upsilon_g

  // TODO: finish summary
*/


// library needed for the multivariate normal distribution
using namespace density;

// parameter transformation, R -> (-1,+1)
template <class Type>
Type bound11(Type x){
	return (1.0-exp(-x))/1.0+exp(-x);
}

template <class Type>
Type invlogitSelectivity(int a, Type mu, Type ups){
	Type tmp = exp((a-(mu-ups/2.0))/ups*10.0);
	return tmp/(1.0+tmp);
} // steepness=10 mimics well selectivitiy curve in (F.7) and (F.8) with ups<5

template <class Type> 
Type square(Type x){
	return x*x;
}

template <class Type>
Type neglogdunif(Type x, Type a, Type b, int n){
	// neg log density of unif(a,b), a<b
	// n = sample size, so that largecst dominates other negloglik contrib
	Type largecst = 10.0*n; // must dominate log(b-a) and other negloglik contrib
	Type res1 = CppAD::CondExpGt(x, a, 0.5*log(b-a), largecst);
	Type res2 = CppAD::CondExpLt(x, b, 0.5*log(b-a), largecst);
	return res1+res2; // neg log density if a<x<b
}



template<class Type>
Type objective_function<Type>::operator() () {

	//----------------------------------------------------------------------------
	// Inputs
	//----------------------------------------------------------------------------

	// Data
	DATA_VECTOR(Ct); // commercial catch, g=1, dim TC
	DATA_VECTOR(S1t); // abundance index survey g=2, dim TS1
	DATA_VECTOR(S2t); // abundance index survey g=3, dim TS2
	DATA_VECTOR(S3t); // abundance index survey g=4, dim TS3

	DATA_MATRIX(patC1); // weighted proportion commercial catch females, dim AxUC
	DATA_MATRIX(patC2); // weighted proportion commercial catch males, dim AxUC
	DATA_VECTOR(ntC); // observed total counts catch, dim UC
	DATA_MATRIX(patS11); // weighted proportion survey 1 females, dim AxUS1
	DATA_MATRIX(patS12); // weighted proportion survey 1 males, dim AxUS1
	DATA_VECTOR(ntS1); // observed total counts survey S1, dim US1

	DATA_IVECTOR(tS1); // C++ indices of S1t values within 0:(TC-1), dim TS1
	DATA_IVECTOR(tS2); // C++ indices of S2t values within 0:(TC-1), dim TS2
	DATA_IVECTOR(tS3); // C++ indices of S3t values within 0:(TC-1), dim TS3
	DATA_IVECTOR(tUC); // C++ indices of patC values within 0:(TC-1), dim UC
	DATA_IVECTOR(tUS1); // C++ indices of patS1 values within 0:(TC-1), dim US1

	DATA_VECTOR(wa1); // average weight females, dim A
	DATA_VECTOR(wa2); // average weight males, dim A
	DATA_VECTOR(ma); // proportion females mature, dim A
	DATA_VECTOR(kappaS1); // known sd of obs error of S1t, dim TS1
	DATA_VECTOR(kappaS2); // known sd of obs error of S2t, dim TS2
	DATA_VECTOR(kappaS3); // known sd of obs error of S3t, dim TS3

	DATA_SCALAR(muS2); // age full selectivity females, S2
	DATA_SCALAR(deltaS2); // age shift full selectivity males, S2
	DATA_SCALAR(upsilonS2); // log scale selectivity curves, S2
	DATA_SCALAR(muS3); // age full selectivity females, S3
	DATA_SCALAR(deltaS3); // age shift full selectivity males, S3
	DATA_SCALAR(upsilonS3); // log scale selectivity curves, S3

	DATA_SCALAR(sigmaR); // sd of proc error of Rt

	DATA_INTEGER(lkhdpropatage); // 1 = Gaussian , 2 = binomial
	DATA_INTEGER(varweight); // 0 = disabled , 1 = inflated var propatage
	DATA_INTEGER(enablepriors); // 0 = disabled , 1 = all enabled


	// Fixed parameters
	PARAMETER(logR0); // log virgin recruitment // TODO: fixed? predicted with Rt?
	PARAMETER(logM1); // log nat mort rate female
	PARAMETER(logM2); // log nat mort rate male
	PARAMETER(logmuC); // log age full selectivity females, catch
	PARAMETER(deltaC); // age shift full selectivity males, catch
	PARAMETER(logupsilonC); // log scale selectivity curves, catch
	PARAMETER(logmuS1); // log age full selectivity females, S1
	PARAMETER(deltaS1); // age shift full selectivity males, S1
	PARAMETER(logupsilonS1); // log scale selectivity curves, S1
	PARAMETER(logh); // log of h in BH SRR
	PARAMETER(logqS1); // log catchabilities S1
	PARAMETER(logqS2); // log catchabilities S2
	PARAMETER(logqS3); // log catchabilities S3
	

	// Random effects
	PARAMETER_VECTOR(logRt); // recruits, dim TC


	//----------------------------------------------------------------------------
	// Setup, procedures and init
	//----------------------------------------------------------------------------

	int TC = Ct.size(); // t = 1, ..., T=74
	int TS1 = S1t.size();
	int TS2 = S2t.size();
	int TS3 = S3t.size();
	int UC = patC1.cols();
	int US1 = patS11.cols();
	int A = patS11.rows(); // a = 1, ..., A=30

	int G = 4; // number of fleets // TODO: config
	int indC = 0; // index for catch data
	int indS1 = 1; // g index for S1
	int indS2 = 2; // g index for S2
	int indS3 = 3; // g index for S3

	Type R0 = exp(logR0);

	vector<Type> Rt = exp(logRt);

	Type M1 = exp(logM1);
	Type M2 = exp(logM2);
	Type expnegM1 = exp(-M1);
	Type expnegM2 = exp(-M2);
	Type sqrtexpnegM1 = sqrt(expnegM1); // exp(-M1*0.5)
	Type sqrtexpnegM2 = sqrt(expnegM2); // exp(-M2*0.5)

	Type muC = exp(logmuC);
	Type upsilonC = exp(logupsilonC);
	Type muS1 = exp(logmuS1);
	Type upsilonS1 = exp(logupsilonS1);

	Type h = exp(logh);

	Type qS1 = exp(logqS1);
	Type qS2 = exp(logqS2);
	Type qS3 = exp(logqS3);

	matrix<Type> Nat1(A,TC); // females
	matrix<Type> Nat2(A,TC); // males

	matrix<Type> sag1(A,G); // selectivity females
	matrix<Type> sag2(A,G); // selectivity males

	vector<Type> Bt(TC); // 1<=t<=T, B_0 is separate

	vector<Type> Vt(TC); // vulnerable biomass
	vector<Type> ut(TC); // alternative to fishing mortality F_t

	matrix<Type> uat1(A,TC); // females
	matrix<Type> uat2(A,TC); // males

	Type nll = 0.0; // init neg loglik


	//----------------------------------------------------------------------------
	// Priors
	//----------------------------------------------------------------------------

	// TODO: config, user specifies hyperparam

	if (enablepriors==1){ // enabled, all values from Table F.4, p. 83
		Type LBR0 = 1.0;
		Type UBR0 = 100000.0;
		nll += neglogdunif(R0, LBR0, UBR0, TC);
		Type meanM1 = 0.07;
		Type sdM1 = 0.007;
		nll -= dnorm(M1, meanM1, sdM1, true);
		nll -= dnorm(M2, meanM1, sdM1, true);
		Type shape1 = 4.573;
		Type shape2 = 2.212;
		nll -= dbeta(h, shape1, shape2, true); // mean=0.674, sd=0.168
		Type LBlogq = -12.0;
		Type UBlogq = 5.0;
		nll += neglogdunif(logqS1, LBlogq, UBlogq, TC);
		nll += neglogdunif(logqS2, LBlogq, UBlogq, TC);
		nll += neglogdunif(logqS3, LBlogq, UBlogq, TC);
		Type meanmuC = 10.5;
		Type sdmuC = 3.15;
		nll -= dnorm(muC, meanmuC, sdmuC, true);
		Type meandeltaC = 0.0;
		Type sddeltaC = 0.3;
		nll -= dnorm(deltaC, meandeltaC, sddeltaC, true);
		Type meanupsilonC = 1.52;
		Type sdupsilonC = 0.456;
		nll -= dnorm(logupsilonC, meanupsilonC, sdupsilonC, true);
		Type meanmuS1 = 13.3;
		Type sdmuS1 = 4.0;
		nll -= dnorm(muS1, meanmuS1, sdmuS1, true);
		Type meandeltaS1 = 0.22;
		Type sddeltaS1 = 0.066;
		nll -= dnorm(deltaS1, meandeltaS1, sddeltaS1, true);
		Type meanupsilonS1 = 2.3; // equivalent to Andy's 3.3
		Type sdupsilonS1 = 1.0;
		nll -= dnorm(logupsilonS1, meanupsilonS1, sdupsilonS1, true);
	}



	//----------------------------------------------------------------------------
	// States dynamics
	//----------------------------------------------------------------------------

	// Selectivities, all fleets
	for (int a = 0; a < A; a++){
			sag1(a,indC) = invlogitSelectivity(a+1, muC, upsilonC); // catch females (F.7)
			sag2(a,indC) = invlogitSelectivity(a+1, muC+deltaC, upsilonC); // catch males (F.8)
			sag1(a,indS1) = invlogitSelectivity(a+1, muS1, upsilonS1); // S1 females (F.7)
			sag2(a,indS1) = invlogitSelectivity(a+1, muS1+deltaS1, upsilonS1); // S1 males (F.8)
			sag1(a,indS2) = invlogitSelectivity(a+1, muS2, upsilonS2); // S2 females (F.7)
			sag2(a,indS2) = invlogitSelectivity(a+1, muS2+deltaS2, upsilonS2); // S2 males (F.8)
			sag1(a,indS3) = invlogitSelectivity(a+1, muS3, upsilonS3); // S3 females (F.7)
			sag2(a,indS3) = invlogitSelectivity(a+1, muS3+deltaS3, upsilonS3); // S3 males (F.8)
	}


	// init Nat, Bt, Rt, Vt, ut and uat
	// Nat, t=1, s=1,2
	for (int a = 0; a < (A-1); a++){
		Nat1(a,0) = 0.5*R0*exp(-M1*a); // (F.4)
		Nat2(a,0) = 0.5*R0*exp(-M2*a); // (F.4)
	}
	Nat1(A-1,0) = 0.5*R0*exp(-M1*(A-1))/(1.0-expnegM1); // (F.5)
	Nat2(A-1,0) = 0.5*R0*exp(-M2*(A-1))/(1.0-expnegM2); // (F.5)
	// Bt, t=0,1
	Type sumwmN1 = 0.0;
	for (int a = 0; a < A; a++){
		sumwmN1 += wa1(a)*ma(a)*Nat1(a,0); // (F.6)
	}
	Bt(0) = sumwmN1; // t=1 (F.6)
	Type B0 = sumwmN1; // t=0 (F.6)
	// Rt, t=1
	Type meanR0 = 4.0*h*R0*B0/((1.0-h)*B0+(5.0*h-1.0)*B0); // (F.10)
	nll -= dnorm(log(Rt(0)), log(meanR0)-square(sigmaR)/2.0, sigmaR, true); // (F.17), (F.18)
	// Vt and ut, t=1
	Type sumwsN1 = 0.0;
	Type sumwsN2 = 0.0;
	for (int a = 0; a < A; a++){
		sumwsN1 += wa1(a)*sag1(a,indC)*Nat1(a,0); // (F.11)
		sumwsN2 += wa2(a)*sag2(a,indC)*Nat2(a,0); // (F.11)
	}
	Vt(0) = sqrtexpnegM1*sumwsN1+sqrtexpnegM2*sumwsN2; // (F.11)
	ut(0) = Ct(0)/Vt(0); // (F.12)
	// uat, t=1
	for (int a = 0; a < A; a++){
		uat1(a,0) = sag1(a,indC)*ut(0); // (F.13)
		uat2(a,0) = sag2(a,indC)*ut(0); // (F.13)
	}


	// dynamics for Rt, Nat, Bt, Vt, ut and uat, 2<=t<=T
	for (int t = 1; t < TC; t++){
		// Rt, based on Bt(t-1)
		Type meanRt = 4.0*h*R0*Bt(t-1)/((1.0-h)*B0+(5.0*h-1.0)*Bt(t-1)); // (F.10)
		nll -= dnorm(log(Rt(t)), log(meanRt)-square(sigmaR)/2.0, sigmaR, true); // (F.17), (F.18)
		// Nat, based on Rt(t) and uat(,t-1)
		Nat1(0,t) = 0.5*Rt(t); // (F.1)
		Nat2(0,t) = 0.5*Rt(t); // (F.1)
		for (int a = 1; a < (A-1); a++){
			Nat1(a,t) = expnegM1*(1.0-uat1(a-1,t-1))*Nat1(a-1,t-1); // (F.2)
			Nat2(a,t) = expnegM2*(1.0-uat2(a-1,t-1))*Nat2(a-1,t-1); // (F.2)
		}
		Nat1(A-1,t) = expnegM1*(1.0-uat1(A-2,t-1))*Nat1(A-2,t-1)
			+ expnegM1*(1.0-uat1(A-1,t-1))*Nat1(A-1,t-1); // (F.3)
		Nat2(A-1,t) = expnegM2*(1.0-uat2(A-2,t-1))*Nat2(A-2,t-1)
			+ expnegM2*(1.0-uat2(A-1,t-1))*Nat2(A-1,t-1); // (F.3)
		// Bt, based on Nat(,t)
		Type sumwmN1 = 0.0;
		for (int a = 0; a < A; a++){
			sumwmN1 += wa1(a)*ma(a)*Nat1(a,t); // (F.9)
		}
		Bt(t) = sumwmN1; // (F.9)
		// Vt and ut, based on Nat(,t)
		Type sumwsN1 = 0.0;
		Type sumwsN2 = 0.0;
		for (int a = 0; a < A; a++){
			sumwsN1 += wa1(a)*sag1(a,indC)*Nat1(a,t); // (F.11)
			sumwsN2 += wa2(a)*sag2(a,indC)*Nat2(a,t); // (F.11)
		}
		Vt(t) = sqrtexpnegM1*sumwsN1+sqrtexpnegM2*sumwsN2; // (F.11)
		ut(t) = Ct(t)/Vt(t); // (F.12)
		// uat, based on ut
		for (int a = 0; a < A; a++){
			uat1(a,t) = sag1(a,indC)*ut(t); // (F.13)
			uat2(a,t) = sag2(a,indC)*ut(t); // (F.13)
		}
	}



	//----------------------------------------------------------------------------
	// Observation equations
	//----------------------------------------------------------------------------

	// survey index S1, t in tS1
	vector<Type> meanS1t(TS1);
	for (int t = 0; t < TS1; t++){
		int ind = tS1(t);
		Type sumuwsN1 = 0.0;
		Type sumuwsN2 = 0.0;
		for (int a = 0; a < A; a++){
			sumuwsN1 += (1.0-0.5*uat1(a,ind))*wa1(a)*sag1(a,indS1)*Nat1(a,ind); // (F.14)
			sumuwsN2 += (1.0-0.5*uat2(a,ind))*wa2(a)*sag2(a,indS1)*Nat2(a,ind); // (F.14)
		}
		meanS1t(t) = qS1*(sqrtexpnegM1*sumuwsN1+sqrtexpnegM2*sumuwsN2); // (F.14)
		nll -= dnorm(log(S1t(t)), log(meanS1t(t))-square(kappaS1(t))/2.0, kappaS1(t), true); // (F.20)
		// nll -= dnorm(log(S1t(t)), log(meanS1t(t)), kappaS1(t), true); // (F.20)
	}

	// survey index S2, t in tS2
	vector<Type> meanS2t(TS2);
	for (int t = 0; t < TS2; t++){
		int ind = tS2(t);
		Type sumuwsN1 = 0.0;
		Type sumuwsN2 = 0.0;
		for (int a = 0; a < A; a++){
			sumuwsN1 += (1.0-0.5*uat1(a,ind))*wa1(a)*sag1(a,indS2)*Nat1(a,ind); // (F.14)
			sumuwsN2 += (1.0-0.5*uat2(a,ind))*wa2(a)*sag2(a,indS2)*Nat2(a,ind); // (F.14)
		}
		meanS2t(t) = qS2*(sqrtexpnegM1*sumuwsN1+sqrtexpnegM2*sumuwsN2); // (F.14)
		nll -= dnorm(log(S2t(t)), log(meanS2t(t))-square(kappaS2(t))/2.0, kappaS2(t), true); // (F.20)
		// nll -= dnorm(log(S2t(t)), log(meanS2t(t)), kappaS2(t), true); // (F.20)
	}

	// survey index S3, t in tS3
	vector<Type> meanS3t(TS3);
	for (int t = 0; t < TS3; t++){
		int ind = tS3(t);
		Type sumuwsN1 = 0.0;
		Type sumuwsN2 = 0.0;
		for (int a = 0; a < A; a++){
			sumuwsN1 += (1.0-0.5*uat1(a,ind))*wa1(a)*sag1(a,indS3)*Nat1(a,ind); // (F.14)
			sumuwsN2 += (1.0-0.5*uat2(a,ind))*wa2(a)*sag2(a,indS3)*Nat2(a,ind); // (F.14)
		}
		meanS3t(t) = qS3*(sqrtexpnegM1*sumuwsN1+sqrtexpnegM2*sumuwsN2); // (F.14)
		nll -= dnorm(log(S3t(t)), log(meanS3t(t))-square(kappaS3(t))/2.0, kappaS3(t), true); // (F.20)
		// nll -= dnorm(log(S3t(t)), log(meanS3t(t)), kappaS3(t), true); // (F.20)
	}

	// prop at age C, t in tUC
	matrix<Type> meanpatC1(A,UC);
	matrix<Type> meanpatC2(A,UC);
	if (lkhdpropatage==1){ // Gaussian lkhd
		Type varcst = 0.0; // ini for correct scope
		if (varweight==0){ // no variance inflation
			varcst = 0.0; // no effect
		} else {
			if (varweight==1){ // inflate prop at age variance by adding varcst
				varcst = 1.0/(10.0*A); // F.5.3, Stanley et al. (2009) // (F.19)
			} else {
				error("varweight must be 0 or 1 (disabled or var inflation).");
			}
		}
		for (int t = 0; t < UC; t++){
			int ind = tUC(t);
			vector<Type> usN1(A); // females
			vector<Type> usN2(A); // males
			Type sumusN1 = 0.0;
			Type sumusN2 = 0.0;
			for (int a = 0; a < A; a++){
				usN1(a) = (1.0-0.5*uat1(a,ind))*sag1(a,indC)*Nat1(a,ind); // (F.15)
				usN2(a) = (1.0-0.5*uat2(a,ind))*sag2(a,indC)*Nat2(a,ind); // (F.15)
				sumusN1 += usN1(a); // (F.15)
				sumusN2 += usN2(a); // (F.15)
			}
			Type sumusN = (sqrtexpnegM1*sumusN1+sqrtexpnegM2*sumusN2); // (F.15)
			for (int a = 0; a < A; a++){
				meanpatC1(a,t) = sqrtexpnegM1*usN1(a)/sumusN; // (F.15)
				meanpatC2(a,t) = sqrtexpnegM2*usN2(a)/sumusN; // (F.15)
				Type sdpat1 = sqrt((patC1(a,t)*(1.0-patC1(a,t))+varcst)/ntC(t)); // (F.19)
				Type sdpat2 = sqrt((patC2(a,t)*(1.0-patC2(a,t))+varcst)/ntC(t)); // (F.19)
				// Type sdpat1 = sqrt((meanpatC1(a,t)*(1.0-meanpatC1(a,t))+varcst)/ntC(t)); // (F.19)
				// Type sdpat2 = sqrt((meanpatC2(a,t)*(1.0-meanpatC2(a,t))+varcst)/ntC(t)); // (F.19)
				nll -= dnorm(patC1(a,t), meanpatC1(a,t), sdpat1, true); // (F.19)
				nll -= dnorm(patC2(a,t), meanpatC2(a,t), sdpat2, true); // (F.19)
			}
		}
	} else {
		if(lkhdpropatage==2){ // binomial lkhd
			for (int t = 0; t < UC; t++){
				int ind = tUC(t);
				vector<Type> usN1(A); // females
				vector<Type> usN2(A); // males
				Type sumusN1 = 0.0;
				Type sumusN2 = 0.0;
				for (int a = 0; a < A; a++){
					usN1(a) = (1.0-0.5*uat1(a,ind))*sag1(a,indC)*Nat1(a,ind); // (F.15)
					usN2(a) = (1.0-0.5*uat2(a,ind))*sag2(a,indC)*Nat2(a,ind); // (F.15)
					sumusN1 += usN1(a); // (F.15)
					sumusN2 += usN2(a); // (F.15)
				}
				Type sumusN = (sqrtexpnegM1*sumusN1+sqrtexpnegM2*sumusN2); // (F.15)
				for (int a = 0; a < A; a++){
					meanpatC1(a,t) = sqrtexpnegM1*usN1(a)/sumusN; // (F.15)
					meanpatC2(a,t) = sqrtexpnegM2*usN2(a)/sumusN; // (F.15)
					nll -= dbinom(patC1(a,t)*ntC(t), ntC(t), meanpatC1(a,t), true);
					nll -= dbinom(patC2(a,t)*ntC(t), ntC(t), meanpatC2(a,t), true);
				}
			}
		} else {
			error("lkhdpropatage must be 1 or 2 (Gaussian or binomial).");
		}
	}

	// prop at age S1, t in tUS1
	matrix<Type> meanpatS11(A,US1);
	matrix<Type> meanpatS12(A,US1);
	if (lkhdpropatage==1){ // Gaussian lkhd
		Type varcst = 0.0; // ini for correct scope
		if (varweight==0){ // no variance inflation
			varcst = 0.0; // no effect
		} else {
			if (varweight==1){ // inflate prop at age variance by adding varcst
				varcst = 1.0/(10.0*A); // F.5.3, Stanley et al. (2009) // (F.19)
			} else {
				error("varweight must be 0 or 1 (disabled or var inflation).");
			}
		}
		for (int t = 0; t < US1; t++){
			int ind = tUS1(t);
			vector<Type> usN1(A); // females
			vector<Type> usN2(A); // males
			Type sumusN1 = 0.0;
			Type sumusN2 = 0.0;
			for (int a = 0; a < A; a++){
				usN1(a) = (1.0-0.5*uat1(a,ind))*sag1(a,indS1)*Nat1(a,ind); // (F.15)
				usN2(a) = (1.0-0.5*uat2(a,ind))*sag2(a,indS1)*Nat2(a,ind); // (F.15)
				sumusN1 += usN1(a); // (F.15)
				sumusN2 += usN2(a); // (F.15)
			}
			Type sumusN = (sqrtexpnegM1*sumusN1+sqrtexpnegM2*sumusN2); // (F.15)
			for (int a = 0; a < A; a++){
				meanpatS11(a,t) = sqrtexpnegM1*usN1(a)/sumusN; // (F.15)
				meanpatS12(a,t) = sqrtexpnegM2*usN2(a)/sumusN; // (F.15)
				Type sdpat1 = sqrt((patS11(a,t)*(1.0-patS11(a,t))+varcst)/ntS1(t)); // (F.19)
				Type sdpat2 = sqrt((patS12(a,t)*(1.0-patS12(a,t))+varcst)/ntS1(t)); // (F.19)
				// Type sdpat1 = sqrt((meanpatS11(a,t)*(1.0-meanpatS11(a,t))+varcst)/ntS1(t)); // (F.19)
				// Type sdpat2 = sqrt((meanpatS12(a,t)*(1.0-meanpatS12(a,t))+varcst)/ntS1(t)); // (F.19)
				nll -= dnorm(patS11(a,t), meanpatS11(a,t), sdpat1, true); // (F.19)
				nll -= dnorm(patS12(a,t), meanpatS12(a,t), sdpat2, true); // (F.19)
			}
		}
	} else {
		if(lkhdpropatage==2){ // binomial lkhd
			for (int t = 0; t < US1; t++){
				int ind = tUS1(t);
				vector<Type> usN1(A); // females
				vector<Type> usN2(A); // males
				Type sumusN1 = 0.0;
				Type sumusN2 = 0.0;
				for (int a = 0; a < A; a++){
					usN1(a) = (1.0-0.5*uat1(a,ind))*sag1(a,indS1)*Nat1(a,ind); // (F.15)
					usN2(a) = (1.0-0.5*uat2(a,ind))*sag2(a,indS1)*Nat2(a,ind); // (F.15)
					sumusN1 += usN1(a); // (F.15)
					sumusN2 += usN2(a); // (F.15)
				}
				Type sumusN = (sqrtexpnegM1*sumusN1+sqrtexpnegM2*sumusN2); // (F.15)
				for (int a = 0; a < A; a++){
					meanpatS11(a,t) = sqrtexpnegM1*usN1(a)/sumusN; // (F.15)
					meanpatS12(a,t) = sqrtexpnegM2*usN2(a)/sumusN; // (F.15)
					nll -= dbinom(patS11(a,t)*ntS1(t), ntS1(t), meanpatS11(a,t), true);
					nll -= dbinom(patS12(a,t)*ntS1(t), ntS1(t), meanpatS12(a,t), true);
				}
			}
		} else {
			error("lkhdpropatage must be 1 or 2 (Gaussian or binomial).");
		}
	}


	//----------------------------------------------------------------------------
	// Outputs
	//----------------------------------------------------------------------------

	ADREPORT(R0);
	ADREPORT(M1);
	ADREPORT(M2);
	ADREPORT(muC);
	ADREPORT(deltaC);
	ADREPORT(upsilonC);
	ADREPORT(muS1);
	ADREPORT(deltaS1);
	ADREPORT(upsilonS1);
	ADREPORT(h);
	ADREPORT(qS1);
	ADREPORT(qS2);
	ADREPORT(qS3);

	// ADREPORT(kappaS1); // fixed, nothing to report
	// ADREPORT(kappaS2); // fixed, nothing to report
	// ADREPORT(kappaS3); // fixed, nothing to report
	// ADREPORT(sigmaR); // fixed, nothing to report

	ADREPORT(Rt); // only strict randeff
	ADREPORT(Bt); // derived quantity
	ADREPORT(Vt); // derived quantity
	ADREPORT(ut); // derived quantity
	ADREPORT(uat1); // derived quantity
	ADREPORT(uat2); // derived quantity
	ADREPORT(sag1); // derived quantity
	ADREPORT(sag2); // derived quantity
	ADREPORT(Nat1); // derived quantity
	ADREPORT(Nat2); // derived quantity

	REPORT(meanS1t);
	REPORT(meanS2t);
	REPORT(meanS3t);
	REPORT(meanpatC1); // mean prop-at-age catch females
	REPORT(meanpatC2); // mean prop-at-age catch males
	REPORT(meanpatS11); // mean prop-at-age S1 females
	REPORT(meanpatS12); // mean prop-at-age S1 males

	return nll;
}

