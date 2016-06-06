#ifndef __ANAL_SPECTRA_CC__
#define __ANAL_SPECTRA_CC__
#include "b3d.h"

//
void CB3D::CalcBalanceQGP(){
	double ETAMAX=parameter::getD(parmap,string("B3D_ETAMAX"),3.0);
	int NRAPIDITYBINS=100;
	double DELY=0.1;
	vector<double> BalNumRapidity_pipi(NRAPIDITYBINS);
	vector<double> BalNumRapidity_pp(NRAPIDITYBINS);
	vector<double> BalNumRapidity_KK(NRAPIDITYBINS);
	vector<double> BalNumRapidity_pK(NRAPIDITYBINS);
	multimap<long int,CPart *> BalancePartMap;
	int ibin,balanceID,balanceMapKey,key,biggestkey=0;
	CPartMap::iterator ppos1,ppos2;
	CPart *part1,*part2;
	double deltarap,rap1,rap2,increment;
	int irap,ID1,ID2;
	parameter::set(parmap,string("B3D_RESONANCES_DECAYS_FILE"),string("progdata/madai/resinfo/decays_pdg_weak.dat"));
	parameter::set(parmap,string("B3D_RESONANCES_INFO_FILE"),string("progdata/madai/resinfo/resonances_pdg_weak.dat"));
	reslist=new CResList();
	reslist->ReadResInfo();
	bool WRITEDETAILS=parameter::getB(parmap,"WRITEDETAILS",false);
	string detailsfilename,dirname,detailsdirname,command;
	string anal_output_filename;
	FILE *anal_output;
	FILE *detailsfile;
	dirname="analysis/"+run_name+"/"+qualifier;
	command="mkdir -p "+dirname;
	system(command.c_str());
	anal_output_filename=dirname+"/results_balanceQGP.dat";
	anal_output=fopen(anal_output_filename.c_str(),"w");
	int nparts,nevents,ID,ipt,iphi,bid1,bid2;
	double pt,mass,pz,rapidity,et,etot,phi;
	double DELYIDITY;
	int neventsmax=parameter::getI(parmap,"B3D_NEVENTSMAX",10);
	CPartMap::iterator ppos,ppos1_upper,ppos1_lower,ppos2_upper,ppos2_lower;
	CPart *part;
	long long int npions=0,nkaons=0,nprotons=0,nomegas=0;
	nevents=0;
	double sigma=0.0;
	int nsigma=0;
	do{
		COLLISIONS=false;
		KillAllParts();
		nparts=ReadBalance(nevents+1);
		printf("nparts=%d\n",nparts);
		printf("PartMap.size=%d DeadPartMap.size=%d FinalPartMap.size=%d\n",int(PartMap.size()),int(DeadPartMap.size()),int(FinalPartMap.size()));
		if(nparts>0){
			nevents+=1;
			PerformAllActions(); // Decays unstable particles
			printf("PartMap.size=%d DeadPartMap.size=%d FinalPartMap.size=%d\n",int(PartMap.size()),int(DeadPartMap.size()),int(FinalPartMap.size()));
			for(ppos=FinalPartMap.begin();ppos!=FinalPartMap.end();ppos++){
				part=ppos->second;
				balanceID=part->balanceID;
				if(balanceID!=-1){
					part=ppos->second;
					pt=sqrt(part->p[1]*part->p[1]+part->p[2]*part->p[2]);
					et=sqrt(part->msquared+pt*pt);
					pz=et*sinh(part->y-part->eta);
					phi=acos(fabs(part->p[1]/pt));
					etot+=sqrt(et*et+part->p[3]*part->p[3]);
					ID=part->resinfo->code;
					if(ID!=111 && abs(ID)!=211 && abs(ID)!=2212 && abs(ID)!=2112 && ID!=22 && abs(ID)!=321 && abs(ID)!=311){
						printf("funny ID=%d\n",ID);
						exit(1);
					}
					balanceMapKey=balanceID;
					if(balanceMapKey>biggestkey)
						biggestkey=balanceMapKey;
					BalancePartMap.insert(CPartPair(balanceMapKey,part));
				}
				else{
					ID=part->resinfo->code;
					//printf("eta=%g, ID=%d\n",part->eta,ID);
					if(abs(ID)==211){
						npions+=1;
					}
					if(abs(ID)==321){
						nkaons+=1;
					}
					if(abs(ID)==2112 || abs(ID)==2212){
						nprotons+=1;
					}
				}
			}
		}
		ppos=BalancePartMap.end();
		ppos--;
		printf("biggestkey=%d\n",int(ppos->first));
		
	}while(!feof(oscarfile) && nevents<neventsmax);
	fclose(oscarfile);
	oscarfile=NULL;
	
	key=0;
	while(key<=biggestkey){
		ppos1_lower=BalancePartMap.lower_bound(key);
		ppos1_upper=BalancePartMap.upper_bound(key);
		ppos2_lower=BalancePartMap.lower_bound(key+1);
		ppos2_upper=BalancePartMap.upper_bound(key+1);
		if(ppos1_lower->first==key && ppos2_lower->first==key+1){
			for(ppos1=ppos1_lower;ppos1!=ppos1_upper;++ppos1){
				part1=ppos1->second;
				bid1=part1->balanceID;
				rap1=part1->y;
				if(part1->balanceID!=key){
					printf("part1 balanceID doesn't match key\n");
					exit(1);
				}
				ID1=part1->resinfo->code;
				//printf("ID1=%d\n",ID1);
				for(ppos2=ppos2_lower;ppos2!=ppos2_upper;++ppos2){
					part2=ppos2->second;
					bid2=part2->balanceID;
					if(bid1!=bid2){
						//printf("ID2=%d\n",ID2);
						if(part2->balanceID!=key+1){
							printf("part2 balanceID doesn't match key\n");
							exit(1);
						}
						ID2=part2->resinfo->code;
						rap2=part2->y;
						deltarap=fabs(rap1-rap2);
						irap=lrint(floor(deltarap/DELY));
						if(irap<NRAPIDITYBINS){
							if(abs(ID1)==211 && abs(ID2)==211){
								if(ID1*ID2<0){
									increment=1.0;
									sigma+=deltarap*deltarap;
									nsigma+=1;
								}
								else
									increment=-1.0;
								//printf("bid1=%d, bid2=%d\n",bid1,bid2);
								BalNumRapidity_pipi[irap]+=increment/DELY;
							}
							if(abs(ID1)==321 && abs(ID2)==321){
								if(ID1*ID2<0)
									increment=1.0;
								else
									increment=-1.0;
								BalNumRapidity_KK[irap]+=increment/DELY;
							}
							if(abs(ID1)==2112 && abs(ID2)==2112){
								if(ID1*ID2<0){
									increment=1.0;
								}
								else
									increment=-1.0;
								BalNumRapidity_pp[irap]+=increment/DELY;
							}
							if(abs(ID1)==2212 && abs(ID2)==2212){
								if(ID1*ID2<0)
									increment=1.0;
								else
									increment=-1.0;
								BalNumRapidity_pp[irap]+=increment/DELY;
							}
						}
					}
				}
			}
		}
		key+=2;
	}
	
	printf("#-------- pipi BF (npions=%lld, nkaons=%lld, nprotons=%lld) ------------\n",
	npions,nkaons,nprotons);
	fprintf(anal_output,"#-------- pipi BF (npions=%lld, nkaons=%lld, nprotons=%lld) ------------\n",
	npions,nkaons,nprotons);
	for(irap=0;irap<NRAPIDITYBINS;irap++){
		printf("%5.2f %g %g %g\n",
		(irap+0.5)*DELY,2*ETAMAX*BalNumRapidity_pipi[irap]/double(npions),2*ETAMAX*BalNumRapidity_KK[irap]/double(nkaons),2*ETAMAX*BalNumRapidity_pp[irap]/double(nprotons));
		fprintf(anal_output,"%5.2f %g %g %g\n",
		(irap+0.5)*DELY,2*ETAMAX*BalNumRapidity_pipi[irap]/double(npions),2*ETAMAX*BalNumRapidity_KK[irap]/double(nkaons),2*ETAMAX*BalNumRapidity_pp[irap]/double(nprotons));
	}
	fclose(anal_output);
	printf("sigma for pions=%g\n",sqrt(sigma/double(nsigma)));
		/*
		if(WRITEDETAILS){
		detailsfilename=detailsdirname+"/balance_pipi.dat";
		detailsfile=fopen(detailsfilename.c_str(),"w");
		fprintf(detailsfile,"# PIPI BALANCE \n");
		fflush(detailsfile);
		fclose(detailsfile);
		}
		*/
		delete reslist;
}

		



#endif
