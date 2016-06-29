#ifndef __ACTION_PERFORM_COLLIDE_CC__
#define __ACTION_PERFORM_COLLIDE_CC__
#include "b3d.h"

void CAction::PerformCollide(){
	int colltype,iproduct,nproducts,wproduct,nrealscatt=0,nfakescatt=0;
	bool productreality=true;
	CPart *part1,*part2,*part,*splitpart;
	CPartMap::iterator ppos;
	CB3DCell *cell;

	ppos=partmap.begin();
	part1=ppos->second;
	++ppos;
	part2=ppos->second;
	part1->actionmother=b3d->nactions;
	part2->actionmother=b3d->nactions;
	
	b3d->GetDeadParts(product);
	colltype=b3d->Collide(part1,part2,nproducts,product,pibsquared);
	if(colltype==0 || nproducts==0){
		b3d->npass+=1;
		return;
	}
	if(colltype==1)
		b3d->nmerge+=1;
	if(colltype==2)
		b3d->nscatter+=1;
	if(colltype==3)
		b3d->ninelastic+=1;
	if(colltype==4)
		b3d->nannihilate+=1;
	
	wproduct=part1->weight*part2->weight;
	if((part1->reality && part2->reality)){
		part1->Kill();
		part2->Kill();
	}
	else{
		productreality=false;
		if(part2->reality){ //make sure part 2 is the fake part and part 1 is the real part
			splitpart=part1;
			part1=part2;
			part2=splitpart;
		}
		nfakescatt=part2->nscatt;
		nrealscatt=part1->nscatt;
		b3d->SplitPart(part1,part2);

		part1->tau_lastint=tau;
		part1->FindActions();

		part2->weight=-part1->weight;
		part2->tau_lastint=tau;
		part2->nscatt=nfakescatt+1;
		part2->FindActions();
	}
		
	for(iproduct=0;iproduct<nproducts;iproduct++){
		part=product[iproduct];
		part->reality=productreality;
		if(productreality)
			part->nscatt=0;
		else
			part->nscatt=nfakescatt+1;
		part->SetMass();
		part->active=true;
		part->weight=wproduct;
		part->tau_lastint=tau;
		part->actionmother=b3d->nactions;
		cell=part->FindCell();
		part->ChangeCell(cell);
		if(part->cell!=NULL){
			if(part->currentmap!=&b3d->PartMap)
				part->ChangeMap(&b3d->PartMap);
		}
		else{
			if(part->currentmap!=&b3d->FinalPartMap)
				part->ChangeMap(&b3d->FinalPartMap);
		}
		part->FindActions();
	}
}

void CAction::PerformCollide_BALANCE(){
	int colltype,iproduct,nproducts,wproduct,nrealscatt=0,nfakescatt=0,balanceID;
	bool productreality=true;
	CPart *part1,*part1a,*part2,*part2a,*part,*splitpart;
	CPartMap::iterator ppos;
	CB3DCell *cell;
	if(part1->weight!=1 || part2->weight!=1){
		printf("part weights don't make sense\n");
		exit(1);
	}

	ppos=partmap.begin();
	part1=ppos->second;
	++ppos;
	part2=ppos->second;
	//make sure part 2 is the fake part and part 1 is the real part
	if(part1->reality && part2->reality){
		printf("Hmmm both particles are real\n");
		printf("balanceIDs = %d, %d\n",part1->balanceID,part2->balanceID);
		exit(1);
	}
	if(!part1->reality && !part2->reality){
		printf("both particles are fake\n");
		exit(1);
	}
	if(part2->reality && !part1->reality){
		splitpart=part1;
		part1=part2;
		part2=splitpart;
	}
	
	part1->actionmother=b3d->nactions;
	part2->actionmother=b3d->nactions;
	b3d->GetDeadParts(product);
	if(!part1->reality || part2->reality){
		printf("wrong reality before Collide\n");
		if(!part1->reality){
			printf("part1 is not real, balanceID=%d\n",part1->balanceID);
			printf("part2 balanceID=%d\n",part2->balanceID);
		}
		if(part2->reality){
			printf("part2 is real, balanceID=%d\n",part2->balanceID);
			printf("part1 balanceID=%d\n",part1->balanceID);
		}
		exit(1);
	}
	//
	colltype=b3d->Collide(part1,part2,nproducts,product,pibsquared);
	//
	if(!part1->reality || part2->reality){
		printf("wrong reality after Collide\n");
		exit(1);
	}
	if(colltype==0 || nproducts==0){
		b3d->npass+=1;
		return;
	}
	
	if(colltype==2){
		//printf("colltype=2, part1->bID=%d, part2->bid=%d\n",part1->balanceID,part2->balanceID);
		part1->tau_lastint=tau;
		part1->nscatt=0;
		part1->balanceID=-1;
		part1->reality=true;
		part1->FindActions();
		if(!part1->reality || part1->balanceID!=-1){
			printf("part1 screwed up\n");
			exit(1);
		}
		part=product[1];
		part->reality=part2->reality;
		part->balanceID=part2->balanceID;
		part->weight=part2->weight;
		part2->Kill();
		part->active=true;
		part->tau_lastint=tau;
		cell=part->FindCell();
		part->ChangeCell(cell);
		if(part->cell!=NULL){
			if(part->currentmap!=&b3d->PartMap)
				part->ChangeMap(&b3d->PartMap);
		}
		else{
			if(part->currentmap!=&b3d->FinalPartMap)
				part->ChangeMap(&b3d->FinalPartMap);
		}
		if(part->reality || part->balanceID==-1){
			printf("shouldn't be real, part->balanceID=%d\n",part->balanceID);
			exit(1);
		}
	}
	else{
		productreality=part2->reality;
		if(productreality){
			printf("productreality should be false\n");
			exit(1);
		}
		balanceID=part2->balanceID;
		nfakescatt=part2->nscatt;
		
		nfakescatt=part2->nscatt;
		b3d->SplitPart(part1,part2);
		part1->tau_lastint=tau;
		part1->FindActions();
		if(part2->resinfo->CheckForNeutral()){
			part2->Kill();
		}
		else{
			part2->resinfo=b3d->reslist->GetResInfoPtr(-part2->resinfo->code);
			part2->weight=1;
			part2->tau_lastint=tau;
			part2->nscatt=nfakescatt+1;
			part2->FindActions();
			part2->reality=false;
			if(part2->balanceID==-1){
				printf("balanceID wrong for part2\n");
				exit(1);
			}
		}
		if(part1->balanceID!=-1 || !part1->reality){
			printf("part1 screwed up\n");
			exit(1);
		}
		
		for(iproduct=0;iproduct<nproducts;iproduct++){
			part=product[iproduct];
			if(!(part->resinfo->CheckForNeutral())){
				part->reality=productreality;
				part->balanceID=balanceID;
				part->nscatt=nfakescatt+1;
				part->SetMass();
				part->active=true;
				part->weight=1;
				part->tau_lastint=tau;
				part->actionmother=b3d->nactions;
				cell=part->FindCell();
				part->ChangeCell(cell);
				if(part->cell!=NULL){
					if(part->currentmap!=&b3d->PartMap)
						part->ChangeMap(&b3d->PartMap);
				}
				else{
					if(part->currentmap!=&b3d->FinalPartMap)
						part->ChangeMap(&b3d->FinalPartMap);
				}
				part->FindActions();
				if(part->reality || part->balanceID==-1){
					printf("products don't make sense\n");
					exit(1);
				}
			}
		}
	}
}

#endif