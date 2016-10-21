#include "sampler.h"
using namespace std;
Csampler* CvolumeElement2D::sampler=NULL;

Csampler::Csampler(CB3D *b3dset){
	b3d=b3dset;
	xyfptr=NULL;
#ifdef __SAMPLER_WRITE_XY__
	xyfptr=fopen("xy.dat","w");
#endif
	VISCOUSCORRECTIONS=true;
	TRIANGLE_FORMAT=true;
	CvolumeElement2D::sampler=this;
	randy=b3d->randy;
	reslist=b3d->reslist;
	cummulative_N=0.0;
	nevents=0;
	cummulative_random=-log(randy->ran());
	Tf=1000.0*parameter::getD(b3d->parmap,"HYDRO_FOTEMP",165.0);
	ETAMAX=parameter::getD(b3d->parmap,"B3D_ETAMAX",1.0);
	NBOSE=parameter::getI(b3d->parmap,"B3D_NBOSE",1);
	TRIANGLE_FORMAT=parameter::getB(b3d->parmap,"SAMPLER_TRIANGLE_FORMAT",true);
	JAKI_FORMAT=parameter::getB(b3d->parmap,"SAMPLER_JAKI_FORMAT",true);
	densityf.clear();
	maxweight.clear();
	boseweight.resize(NBOSE+1);
	reslist->CalcEoS(Tf,epsilonf,Pf,nhadronsf,densityf,boseweight,maxweight);
	sf=(epsilonf+Pf)/Tf;
	lambdaf=GetLambda(Tf,Pf,epsilonf);
}

Csampler::~Csampler(){
	if(xyfptr!=NULL){
		fclose(xyfptr);
		xyfptr=NULL;
	}
}

int Csampler::MakeB3DEvent(){
	double Omega0Sum=0.0;
	nevents+=1;
	int ielement,nparts=0;
	for(ielement=0;ielement<nelements;ielement++){
		Omega0Sum+=element[ielement].Omega[0];
		nparts+=element[ielement].MakeParts();
	}
	printf("Event %4d sampling: %d parts created\n",nevents,int(b3d->PartMap.size()));
	return nparts;
}

int Csampler::MakeB3DEvent_UniformXY(){
	CvolumeElement2D elem;
	elem.epsilon=epsilonf;
	elem.P=Pf;
	elem.nhadrons=nhadronsf;
	elem.density=&densityf;
	elem.T=Tf;
	elem.lambda=lambdaf;
	nevents+=1;
	nparts=elem.MakeParts_UniformXY();
	printf("Event %4d sampling: %d parts created\n",nevents,int(b3d->PartMap.size()));
	return nparts;
}

void Csampler::CalcPiFromParts(){
	double pi[4][4]={0.0};
	CPartMap::iterator ppos;
	CPart *part;
	double *p,pressure,Ptest=0.0,etest=0.0;
	int ipart,alpha,beta;
	double volume=nelements*element[0].Omega[0];
	for(ppos=b3d->DeadPartMap.begin();ppos!=b3d->PartMap.end();ppos++){
		part=ppos->second;
		pressure=(part->p[1]*part->p[1]+part->p[2]*part->p[2]+part->p[3]*part->p[3])/(3.0*part->p[0]);
		Ptest+=pressure;
		etest+=p[0];
		for(alpha=1;alpha<4;alpha++){
			pi[alpha][alpha]-=pressure;
			for(beta=1;beta<4;beta++){
				pi[alpha][beta]+=p[alpha]*p[beta]/p[0];
			}
		}
	}
	printf("----From particle momenta  -----\n");
	printf("epsilon=%10.4f, P=%10.4f\n",etest/volume,Ptest/volume);
	printf("pi_ij/P\n");
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++){
			pi[alpha][beta]=pi[alpha][beta]/volume;
			printf("%9.6f ",pi[alpha][beta]/element[0].P);
		}
		printf("\n");
	}
	printf("---- From pi_ij  input -------\n");
	printf("epsilon=%10.4f, P=%10.4f\n",element[0].epsilon,element[0].P);
	printf("pi_ij/P\n");
	for(alpha=1;alpha<4;alpha++){
		for(beta=1;beta<4;beta++){
			printf("%9.6f ",element[0].pitilde[alpha][beta]/element[0].P);
		}
		printf("\n");
	}
}

void Csampler::ReadVolumeElements2D(){
	if(TRIANGLE_FORMAT)
		ReadVolumeElements2D_triangles();
	else if(JAKI_FORMAT)
		ReadVolumeElements2D_Jaki();
	else{
		printf("Sampler format screwed up\n");
		exit(1);
	}
}
	
void Csampler::ReadVolumeElements2D_triangles(){
	string filename;
	CvolumeElement2D *elem;
	double tau,deltau,H;
	double pixx,pixy,pixz,piyy,piyz,pizz;
	double dTxxoverH,dTxyoverH,dTxzoverH,dTyyoverH,dTyzoverH,dTzzoverH;
	int ivertex,iv1,iv2,iv3,ielement,initarraysize=1000;
	char dummy[300];
	
	vertex.clear();
	filename="model_output/"+b3d->run_name+"/"+b3d->qualifier+"/vertices2D.dat";
	FILE *fptr=fopen(filename.c_str(),"r");
	
	fscanf(fptr,"%d",&ivertex);
	while(!feof(fptr)){
		if(ivertex==vertex.size())
			vertex.resize(vertex.size()+initarraysize);
		fscanf(fptr,"%lf %lf %lf %lf %lf",&vertex[ivertex].r[0],&vertex[ivertex].r[1],&vertex[ivertex].r[2],
		&vertex[ivertex].ux,&vertex[ivertex].uy);
		fscanf(fptr,"%d",&ivertex);
	}
	fclose(fptr);
	nvertices=vertex.size();

	element.clear();
	nelements=0;
	filename="model_output/"+b3d->run_name+"/"+b3d->qualifier+"/triangles2D.dat";
	fptr=fopen(filename.c_str(),"r");
	
	ielement=0;
	double e0,e1,e2;
	fscanf(fptr,"%d %d %d",&iv1,&iv2,&iv3);
	double pi33overPbar=0.0,pi33overPbarnorm=0.0;
	while(!feof(fptr)){
		if(element.size()==ielement)
			element.resize(element.size()+initarraysize);
		elem=&element[ielement];
		fscanf(fptr,"%lf %lf %lf",
		&(elem->Omega[0]),&(elem->Omega[1]),&(elem->Omega[2]));

		fscanf(fptr,"%lf %lf",&elem->T,&elem->Xscale);
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf",
		&dTxxoverH,&dTxyoverH,&dTxzoverH,&dTyyoverH,&dTyzoverH,&dTzzoverH);
		elem->vertex[0]=&vertex[iv1];
		elem->vertex[1]=&vertex[iv2];
		elem->vertex[2]=&vertex[iv3];
		elem->epsilon=epsilonf;
		elem->density=&densityf;
		elem->P=Pf;
		elem->lambda=lambdaf;
		elem->nhadrons=nhadronsf;
		H=(epsilonf+Pf);
		//elem->FillOutShearTensor(H*dTxxoverH,H*dTxyoverH,H*dTxzoverH,H*dTyyoverH,H*dTyzoverH,H*dTzzoverH);
		elem->pitilde[1][1]=H*dTxxoverH;
		elem->pitilde[1][2]=elem->pitilde[2][1]=H*dTxyoverH;
		elem->pitilde[1][3]=elem->pitilde[3][1]=H*dTxzoverH;
		elem->pitilde[2][2]=H*dTyyoverH;
		elem->pitilde[2][3]=elem->pitilde[3][2]=H*dTyzoverH;
		elem->pitilde[3][3]=H*dTzzoverH;
		elem->CalcOmegamax();
		pi33overPbar+=H*dTzzoverH*nhadronsf*elem->Omegamax/Pf;
		pi33overPbarnorm+=nhadronsf*elem->Omegamax;
		ielement+=1;
		fscanf(fptr,"%d %d %d",&iv1,&iv2,&iv3);
	}
	pi33overPbar=pi33overPbar/pi33overPbarnorm;
	nelements=ielement;
	printf("data read in, nelements=%d, pi33/P=%g\n",nelements,pi33overPbar);
}

void Csampler::ReadVolumeElements2D_Jaki(){
	string filename;
	CvolumeElement2D *elem;
	double dumbo,udotn,PIbulk;
	int alpha,beta;
	double sigma,PI,**pi,*u;
	u=new double[4];
	pi=new double *[4];
	for(alpha=0;alpha<4;alpha++)
		pi[alpha]=new double[4];
	int ielement,initarraysize=1000;
	char dummy[300];
	element.clear();
	nelements=0;
	filename="model_output/"+b3d->run_name+"/"+b3d->qualifier+"/hydro_Jaki2D.dat";
	printf("opening %s\n",filename.c_str());
	FILE *fptr=fopen(filename.c_str(),"r");
	ielement=0;

	while(!feof(fptr)){
		if(element.size()==ielement)
			element.resize(element.size()+initarraysize);
		elem=&element[ielement];
		fscanf(fptr,"%lf %lf %lf",&(elem->Omega[0]),&(elem->Omega[1]),&(elem->Omega[2]));
		fscanf(fptr,"%lf %lf %lf %lf",&dumbo,&(elem->ux),&(elem->uy),&sigma);
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf",
		&PIbulk,&pi[0][0],&pi[1][1],&pi[2][2],&pi[3][3],&pi[1][2]);
		fscanf(fptr,"%lf %lf %lf",&(elem->tau),&(elem->x),&(elem->y));
		elem->Xscale=1.0; elem->Omega[3]=0.0; elem->T=155.0; // check Temperature!
		if(!feof(fptr)){
			u[1]=elem->ux; u[2]=elem->uy;
			u[0]=sqrt(1.0+u[1]*u[1]+u[2]*u[2]); u[3]=0.0;
			udotn=u[0]*elem->Omega[0]-u[1]*elem->Omega[1]-u[2]*elem->Omega[2];
			for(alpha=0;alpha<3;alpha++)
				elem->Omega[alpha]=elem->Omega[alpha]*sigma/udotn;
		
			pi[3][3]=pi[3][3]*(elem->tau)*(elem->tau);
		
			pi[2][1]=pi[1][2];
			pi[2][3]=pi[3][2]=pi[1][3]=pi[3][1]=pi[0][3]=pi[3][0]=0.0;
			pi[0][1]=(pi[1][1]*u[1]+pi[1][2]*u[2])/u[0];
			pi[1][0]=pi[0][1];
			pi[0][2]=(pi[2][1]*u[1]+pi[2][2]*u[2])/u[0];
			pi[2][0]=pi[0][2];
			pi[0][3]=pi[3][0]=0.0;
			pi[3][3]=-pi[1][1]-pi[2][2]+pi[0][0]; 
			pi[0][0]=(pi[0][1]*u[1]+pi[0][2]*u[2])/u[0];

			for(alpha=0;alpha<4;alpha++){
				for(beta=0;beta<4;beta++)
					pi[alpha][beta]*=HBARC;
			}
		
			elem->epsilon=epsilonf;
			elem->density=&densityf;
			elem->P=Pf;
			elem->lambda=lambdaf;
			elem->nhadrons=nhadronsf;
			elem->CalcOmegamax();
			ielement+=1;
			//printf("Omega=(%g,%g,%g,%g),sigma=%g, Omegamax=%g\n",
			//elem->Omega[0],elem->Omega[1],elem->Omega[2],elem->Omega[3],sigma,elem->Omegamax);
			//Misc::Pause();
		}
	}
	nelements=ielement;
	printf("Exiting happily\n");
}

void Csampler::ReadVolumeElements3D(){
	string filename;
	CvolumeElement2D *elem;
	double pixx, pixy, pixz, piyy, piyz, pizz;
	
	int ielement, initarraysize=1000;
	element.clear();
	ielement=0;
	nelements=0;

	filename="model_output/"+b3d->run_name+"/"+b3d->qualifier+"/hydro3D.dat";
	printf("\n***opening %s***\n\n", filename.c_str());
	FILE *fptr=fopen(filename.c_str(),"r");

	while(!feof(fptr)){
		if(element.size()==ielement)
			element.resize(element.size()+initarraysize);
		elem=&element[ielement];

		fscanf(fptr,"%lf %lf %lf %lf",&elem->tau,&elem->eta,&elem->x,&elem->y);
		fscanf(fptr,"%lf %lf %lf",&(elem->ux),&(elem->uy),&(elem->uz));
		fscanf(fptr,"%lf %lf %lf %lf",&(elem->Omega[0]),&(elem->Omega[1]),&(elem->Omega[2]),&(elem->Omega[3]));
		fscanf(fptr,"%lf %lf %lf %lf %lf %lf",&pixx,&pixy,&pixz,&piyy,&piyz,&pizz);
		
		printf("%lf %lf %lf %lf",elem->tau,elem->eta,elem->x,elem->y);
		printf("%lf %lf %lf",elem->ux,elem->uy,elem->uz);
		printf("%lf %lf %lf %lf",elem->Omega[0],elem->Omega[1],elem->Omega[2],elem->Omega[3]);
		printf("%lf %lf %lf %lf %lf %lf\n",pixx,pixy,pixz,piyy,piyz,pizz);
		
		elem->epsilon=epsilonf;
		elem->density=&densityf;
		elem->P=Pf;
		elem->lambda=lambdaf;
		elem->nhadrons=nhadronsf;

		elem->pitilde[1][1]=pixx;
		elem->pitilde[1][2]=elem->pitilde[2][1]=pixy;
		elem->pitilde[1][3]=elem->pitilde[3][1]=pixz;
		elem->pitilde[2][2]=piyy;
		elem->pitilde[2][3]=elem->pitilde[3][2]=piyz;
		elem->pitilde[3][3]=pizz;

		printf("[%lf %lf %lf]\n",elem->pitilde[1][1],elem->pitilde[1][2],elem->pitilde[1][3]);
		printf("[%lf %lf %lf]\n",elem->pitilde[2][1],elem->pitilde[2][2],elem->pitilde[2][3]);
		printf("[%lf %lf %lf]\n",elem->pitilde[3][1],elem->pitilde[3][2],elem->pitilde[3][3]);

		printf("ielement=%d\n",ielement);
		ielement+=1;
	}
	nelements=ielement;
}

double Csampler::GetLambda(double T,double P,double epsilon){
	int iQ,n,i;
	const int nmax=70;
	double G[nmax+5];
	double m,degen,z,Ipp=0.0,Ipptest=0.0,dIpp,Ptest=0.0,J,nfact,sign,alpha;
	double dIpptest=0.0,dp=4.0,p,e,lambdafact;
	CResInfo *resinfo;
	CResInfoMap::iterator rpos;
	for(rpos=reslist->resmap.begin();rpos!=reslist->resmap.end();rpos++){
		resinfo=rpos->second;
		if(resinfo->code!=22){
			m=resinfo->mass;
			degen=(2.0*resinfo->spin+1);
			z=m/T;
			alpha=0.0;

			G[0]=gsl_sf_gamma_inc(5,z)*pow(z,-5);
			for(int i=1;i<nmax+5;i++){
				n=5-2*i;
				if(n!=-1)	G[i]=(-exp(-z)/n)+(G[i-1]*z*z-z*exp(-z))/((n+1.0)*n);
				else G[i]=gsl_sf_gamma_inc(-1,z)*z;
			}
			J=0.0;
			nfact=1.0;
			sign=1.0;
			for(n=0;n<nmax;n+=1){
				if(n>0) sign=-1.0;
				J+=sign*nfact*(G[n]-2.0*G[n+1]+G[n+2]);
				nfact=nfact*0.5/(n+1.0);
				if(n>0) nfact*=(2.0*n-1.0);
			}
			dIpp=degen*exp(alpha)*pow(m,4)*(-z*J+15.0*gsl_sf_bessel_Kn(2,z)/(z*z));
			dIpp=dIpp/(60.0*PI*PI*HBARC*HBARC*HBARC);
			/*
			dIpptest=0.0;
			for(p=0.5*dp;p<3000;p+=dp){
			e=sqrt(m*m+p*p);
			dIpptest+=degen*(4.0*PI/(8.0*PI*PI*PI*HBARC*HBARC*HBARC))*p*p*dp*exp(-e/T)*( (2.0/3.0)*(p*p/e) - (2.0/15.0)*pow(p,4)/pow(e,3) );
			//dIpptest+=degen*(4.0*PI/(8.0*PI*PI*PI*HBARC*HBARC*HBARC))*p*p*dp*exp(-e/T)*((2.0/15.0)*pow(p,4)/pow(e,3) );
			}
			Ipptest+=dIpptest;
			*/
			
			Ipp+=dIpp;
			//Ptest+=Ipptest;
			//printf("dIpptest=%g =? %g, ratio=%g\n",dIpptest,dIpp,dIpptest/dIpp);
		}
	}
	lambdafact=(2.0*P-4.0*Ipp)/(P+epsilon);
	//printf("P=%g, epsilon=%g\n",P,epsilon);
	//printf("Ipp=%g =? %g, lambdafact=%g\n",Ipptest,2.0*P-4.0*Ipp,lambdafact);
	return lambdafact;
}


