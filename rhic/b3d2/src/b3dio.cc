#ifndef __B3DIO_CC__
#define __B3DIO_CC__
#include "b3d.h"

double CB3D::WriteOSCAR(int ievent){
	CB3DBinaryPartInfo bpart;
	char dummy[100];
	if(oscarfile==NULL){
		if(BINARY_RW)
			oscarfile=fopen(oscarfilename.c_str(),"wb");
		else{
			oscarfile=fopen(oscarfilename.c_str(),"w");
			fprintf(oscarfile,":OSCAR1997a\n");
			fprintf(oscarfile,"ipart -- id -- p[4] -- m -- x[4]\n");
			fprintf(oscarfile,"b3d output\n");
		}
	}
	int nparts=PartMap.size()+FinalPartMap.size();
	if(BINARY_RW){
		fwrite(&ievent,sizeof(int),1,oscarfile);
		fwrite(&nparts,sizeof(int),1,oscarfile);
	}
	else
		fprintf(oscarfile,"%7d %6d    %8.5f     %8.5f\n",ievent,nparts,parameter::getD(parmap,"GLAUBER_B",0.0),
	parameter::getD(parmap,"GLAUBER_B",0.0));
	double v,pperp,eperp,dnchdeta=0.0,dnchdy=0,t,twrite,tauwrite,etawrite,eta,deleta,y;
	FourVector rwrite,pwrite;
	double mass;
	int ipart,nmesons=0,nch=0;
	CPart *part;
	CPartMap::iterator ppos;
	ipart=0;
	ppos=FinalPartMap.begin();
	do{
		if(ppos==FinalPartMap.end())
			ppos=PartMap.begin();
		part=ppos->second;
		part->Propagate(part->tau_lastint);
		if(BJORKEN && fabs(part->eta)>ETAMAX){
			part->CyclicReset();
		}
		if(part->resinfo->baryon!=0)
			nbaryons+=1;
		if(BINARY_RW){
			bpart.ID=part->resinfo->code;
			bpart.tau=part->tau0;
			bpart.x=part->r[1];
			bpart.y=part->r[2];
			bpart.eta=part->eta;
			bpart.px=part->p[1];
			bpart.py=part->p[2];
			bpart.rapidity=part->y;
			bpart.weight=part->weight;
			bpart.reality=part->reality;
			fwrite(&bpart,sizeof(bpart),1,oscarfile);
		}
		else
			fprintf(oscarfile,"%5d %5d %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %12.6e %d %d\n",
		ipart,part->resinfo->code,part->p[1],part->p[2],part->p[3],part->p[0],sqrt(part->msquared),part->r[1],part->r[2],part->r[3],part->r[0],part->weight,int(part->reality));
		if(ppos==PartMap.end()){
			printf("ppos shouldn't be here\n");
		}
		++ppos;
		ipart+=1;
	} while(ipart<nparts);
	return dnchdy/(2.0*ETAMAX);
}

void CB3D::ReadOSCARHeader(){
	int ndead=3,idead;
	char dummy[200];
	oscarfilename="model_output/"+run_name+"/"+qualifier+"/oscar.dat";
	if(BINARY_RW)
		oscarfile=fopen(oscarfilename.c_str(),"rb");
	else{
		oscarfile=fopen(oscarfilename.c_str(),"r");
		for(idead=0;idead<ndead;idead++)
			fgets(dummy,200,oscarfile);
	}
}

int CB3D::ReadOSCAR(int ievent){
	CB3DBinaryPartInfo bpart;
	CResInfo *resinfo;
	double p[4],r[4],mass,rapidity,eta,tau0;
	int weight,reality_int,ID;
	int nparts,nparts_read,ipart=0;
	int ievent_read;
	bool reality;
	double bmin,bmax; // impact parameter
	double mtot,mothermass;
	CPart *mother;
	tau=0.0;
	if(oscarfile==NULL){
		ReadOSCARHeader();
	}
	if(BINARY_RW){
		fread(&ievent_read,sizeof(int),1,oscarfile);
		fread(&nparts_read,sizeof(int),1,oscarfile);
	}
	else{
		fscanf(oscarfile,"%d %d %lf %lf",&ievent_read,&nparts_read,&bmin,&bmax);
		if(!feof(oscarfile) && ievent_read!=ievent){
			printf("trying to read wrong event, ievent=%d, ievent_read=%d\n",ievent,ievent_read);
			exit(1);
		}
	}
	if(feof(oscarfile)){
		return 0;
	}
	for(ipart=0;ipart<nparts_read;ipart++){
		mother=GetDeadPart();
		if(BINARY_RW){
			fread(&bpart,sizeof(bpart),1,oscarfile);
			ID=bpart.ID;
			tau0=bpart.tau;
			r[1]=bpart.x;
			r[2]=bpart.y;
			eta=bpart.eta;
			p[1]=bpart.px;
			p[2]=bpart.py;
			rapidity=bpart.rapidity;
			resinfo=reslist->GetResInfoPtr(ID);
			mass=resinfo->mass;
			weight=bpart.weight;
			reality=bpart.reality;
		}
		else{
			fscanf(oscarfile,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d",
			&ipart,&ID,&p[1],&p[2],&p[3],&p[0],&mass,&r[1],&r[2],&r[3],&r[0],&weight,&reality_int);
			tau0=sqrt(r[0]*r[0]-r[3]*r[3]);
			eta=asinh(r[3]/tau0);
			rapidity=asinh(p[3]/p[0]);
			reality=true;
			if(reality_int==0)
				reality=false;
		}
		mother->Init(ID,r[1],r[2],tau0,eta,p[1],p[2],mass,rapidity,weight,reality);
	}
	return nparts_read;
}

void CB3D::ReadBalanceParts(){
	vector<int> netcharge(50000,0);
	string filename=parameter::getS(parmap,"BALANCE_INPUT_FILENAME","vinzentdata/hadrons.csv");
	printf("reading %s\n",filename.c_str());
	FILE *fptr=fopen(filename.c_str(),"r");
	int nqgppions=0;
	double sigma=0;
	int nsigma=0;
	CPartMap::iterator ppos;
	vector<CPart *> newpart;
	newpart.reserve(10);
	CResInfo *resinfo;
	char dummy[120];
	fgets(dummy,120,fptr);
	double x,y,z,t,E,px,py,pz,mt,tau,eta,eta0,rapidity,mass;
	int ibalance,pid,intweight=1,ibalread,iibalread,ibalpair,oldibalpair=-1;
	int ibp,nbalpairs=0,nparts=0;
	bool reality=false,paircheck,evencheck,oddcheck;
	CPart *part;
	ibalance=0;
	paircheck=false;
	do{
		fscanf(fptr,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf",&pid,&ibalread,&x,&y,&z,&tau,&E,&px,&py,&pz);
		eta0=asinh(z/t);
		resinfo=reslist->GetResInfoPtr(pid);
		if(ibalread==netcharge.size())
			netcharge.resize(netcharge.size()+5000);
		for(iibalread=iibalread;iibalread<netcharge.size();iibalread++)
			netcharge[iibalread]=0;
		netcharge[ibalread]+=resinfo->charge;
		mass=resinfo->mass;
		mt=sqrt(mass*mass+px*px+py*py);
		E=sqrt(mt*mt+pz*pz);
		eta=1.0*randy->gauss();
		if(abs(pid)==211){
			nsigma+=1;
			sigma+=eta*eta;
		}
		z=tau*sinh(eta);
		t=tau*cosh(eta);
		rapidity=(eta-eta0)+0.5*log((E+pz)/(E-pz));
		E=mt*cosh(rapidity);
		pz=mt*sinh(rapidity);
		//t=sqrt(tau*tau+z*z);
		//eta=atanh(z/t);
		if(abs(pid)==211 && ibalread!=-1)
			nqgppions+=1;
		if(!feof(fptr)){
			ibalpair=lrint(floor(double(ibalread)/2.0));
			if(ibalpair!=oldibalpair){
				//if previous round of ibalpair didn't have any good pairs, kill all particles
				if(paircheck==false){
					for(nparts=0;nparts<newpart.size();nparts++){
						newpart[nparts]->Kill();
					}
					newpart.clear();
					nparts=0;
				}
				newpart.clear();
				//printf("--------------------\n");
				nparts=0;
				paircheck=false;
				evencheck=false;
				oddcheck=false;
			}
			if(ibalread%2==0)
				evencheck=true;
			if(ibalread%2==1)
				oddcheck=true;
			if(evencheck && oddcheck){
				if(!paircheck)
					nbalpairs+=1;
				paircheck=true;
			}
			//printf("ibalread=%d, eta=%g\n",ibalread,eta);
			newpart.push_back(GetDeadPart());
			rapidity=0.5*log((E+pz)/(E-pz));
			//tau=sqrt(t*t-z*z);
			//eta=asinh(z/tau);
			if(t<z){
				printf("t<z???, tau=%g\n",tau);
				exit(1);
			}
			resinfo=reslist->GetResInfoPtr(pid);
			mass=resinfo->mass;
			if(ibalread!=-1)
				reality=false;
			else
				reality=true;
			newpart[nparts]->InitBalance(pid,x,y,tau,eta,px,py,mass,rapidity,intweight,reality,ibalread);
			if(newpart[nparts]->resinfo->CheckForNeutral()){
				printf("why are we generating a neutral particle?\n");
				exit(1);
			}
			nparts+=1;
		}
		oldibalpair=ibalpair;
	}while(!feof(fptr));
	int netbal=0;
	for(ibalread=0;ibalread<netcharge.size();ibalread+=2){
		netbal+=netcharge[ibalread]*netcharge[ibalread+1];
	}
	printf("netbal=%d\n",netbal);
	printf("--- nbalpairs=%d, nqgppions=%d, sigma=%g\n",nbalpairs,nqgppions,sqrt(sigma/double(nsigma)));
}

double CB3D::WriteBalance(int ievent){
	printf("In WriteBalance..... FinalPartMap.size()=%d, PartMap.size()=%d\n",int(FinalPartMap.size()),int(PartMap.size()));
	CB3DBinaryBalancePartInfo bpart;
	double sigma=0;
	int nsigma=0;
	char dummy[100];
	if(oscarfile==NULL){
		if(BINARY_RW)
			oscarfile=fopen(oscarfilename.c_str(),"wb");
		else{
			oscarfile=fopen(oscarfilename.c_str(),"w");
			fprintf(oscarfile,":OSCAR1997a\n");
			fprintf(oscarfile,"ipart -- id -- p[4] -- m -- x[4]\n");
			fprintf(oscarfile,"b3d output\n");
		}
	}
	int nparts=PartMap.size()+FinalPartMap.size();
	if(BINARY_RW){
		fwrite(&ievent,sizeof(int),1,oscarfile);
		fwrite(&nparts,sizeof(int),1,oscarfile);
	}
	else
		fprintf(oscarfile,"%7d %6d    %8.5f     %8.5f\n",ievent,nparts,parameter::getD(parmap,"GLAUBER_B",0.0),
	parameter::getD(parmap,"GLAUBER_B",0.0));
	double v,pperp,eperp,dnchdeta=0.0,dnchdy=0,t,twrite,tauwrite,etawrite,eta,deleta,y;
	FourVector rwrite,pwrite;
	double mass,rapidity;
	int ipart,nmesons=0,nch=0;
	CPart *part;
	CPartMap::iterator ppos;
	ipart=0;
	ppos=FinalPartMap.begin();
	do{
		if(ppos==FinalPartMap.end())
			ppos=PartMap.begin();
		part=ppos->second;
		part->Propagate(part->tau_lastint);
		if(BJORKEN && fabs(part->eta)>ETAMAX){
			part->CyclicReset();
		}
		if(part->resinfo->baryon!=0)
			nbaryons+=1;
		if(BINARY_RW){
			bpart.ID=part->resinfo->code;
			bpart.balanceID=part->balanceID;
			bpart.px=part->p[1];
			bpart.py=part->p[2];
			bpart.rapidity=part->y;
			sigma+=part->y*part->y;
			nsigma+=1;
			//if(bpart.balanceID!=-1)
			//printf("final part in WriteBalance, eta=%g\n",eta);
			//if(bpart.balanceID!=-1)
			fwrite(&bpart,sizeof(bpart),1,oscarfile);
		}
		else{
			rapidity=0.5*log((part->p[0]+part->p[3])/(part->p[0]-part->p[3]));
			if(bpart.balanceID!=-1){
				fprintf(oscarfile,"%5d %5d %7d %12.6e %12.6e %12.6e\n",
				ipart,part->resinfo->code,part->balanceID,part->p[1],part->p[2],rapidity);
			}
		}
		++ppos;
		ipart+=1;
	} while(ipart<nparts);
	printf("WriteBalance -- sigma=%g\n",sqrt(sigma/double(nsigma)));
	return dnchdy/(2.0*ETAMAX);
}

int CB3D::ReadBalance(int ievent){
	CB3DBinaryBalancePartInfo bpart;
	CResInfo *resinfo;
	double p[4],r[4],mass,rapidity,eta,tau0;
	int weight,reality_int,ID,balanceID;
	int nparts,nparts_read,ipart=0;
	int ievent_read;
	bool reality;
	double sigma=0;
	int nsigma=0;
	double bmin,bmax; // impact parameter
	double mtot,mothermass;
	CPart *mother;
	tau=0.0;
	if(oscarfile==NULL){
		ReadOSCARHeader();
	}
	if(BINARY_RW){
		fread(&ievent_read,sizeof(int),1,oscarfile);
		fread(&nparts_read,sizeof(int),1,oscarfile);
	}
	else{
		fscanf(oscarfile,"%d %d %lf %lf",&ievent_read,&nparts_read,&bmin,&bmax);
		if(!feof(oscarfile) && ievent_read!=ievent){
			printf("trying to read wrong event, ievent=%d, ievent_read=%d\n",ievent,ievent_read);
			exit(1);
		}
	}
	if(feof(oscarfile)){
		return 0;
	}
	for(ipart=0;ipart<nparts_read;ipart++){
		mother=GetDeadPart();
		if(BINARY_RW){
			fread(&bpart,sizeof(bpart),1,oscarfile);
			ID=bpart.ID;
			p[1]=bpart.px;
			p[2]=bpart.py;
			rapidity=bpart.rapidity;
			balanceID=bpart.balanceID;
		}
		else{
			printf("should only work for binary rw\n");
			exit(1);
			//fscanf(oscarfile,"%d %d %lf %lf %lf",&ipart,&ID,&p[1],&p[2],&rapidity);
		}
		//printf("read in ipart=%d, ID=%d, p=(%g,%g), y=%g, nparts_read=%d\n",
		//ipart,ID,p[1],p[2],rapidity,nparts_read);
		r[1]=r[2]=0.0;
		eta=rapidity;
		tau0=10.0;
		weight=1.0;
		reality=true;
		resinfo=reslist->GetResInfoPtr(ID);
		mass=resinfo->mass;
		if(resinfo->charge!=0 || resinfo->decay){
			mother->InitBalance(ID,r[1],r[2],tau0,eta,p[1],p[2],mass,rapidity,weight,reality,balanceID);
			sigma+=rapidity*rapidity;
			nsigma+=1;
		}
	}
	printf("ReadBalance -- sigma=%g\n",sqrt(sigma/double(nsigma)));
	return nparts_read;
}

void CB3D::WriteDens(){
	string densfilename="model_output/"+run_name+"/"+qualifier+"/dens.dat";
	FILE *densfile = fopen(densfilename.c_str(),"w");
	fprintf(densfile,"#ix iy  dens[itau=0] dens[itau=1]...\n");
	double dxy;
	int ix,iy,ieta,itau;
	for(ix=0;ix<2*NXY;ix++){
		for(iy=0;iy<2*NXY;iy++){
			fprintf(densfile,"%3d %3d",ix,iy);
			for(itau=0;itau<DENSWRITE_NTAU;itau++){
				dxy=0.0;
				for(ieta=0;ieta<2*NETA;ieta++){
					dxy+=cell[ix][iy][ieta]->dens[itau];
				}
				fprintf(densfile," %6.0f",dxy);
			}
			fprintf(densfile,"\n");
		}
	}
	fclose(densfile);
}

#endif
