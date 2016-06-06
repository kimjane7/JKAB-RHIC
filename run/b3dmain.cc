#include "b3d.h"
#include "qualifier.h"

using namespace std;

int main(int argc, char *argv[]){
	if (argc != 2) {
		printf("Usage: b3d run_name\n");
		exit(-1);
  }
	double dnchdy=0.0,tau;
	int nparts;
	long long int npartstot;
	long long int ncolls=0,nannihilate=0,nregen=0,nbaryons=0,norm;
	int ievent=0,iqual,nevents;
	string run_name=argv[1];
	CB3D *b3d=new CB3D(run_name);
	b3d->InitCascade();
	CQualifiers qualifiers;
	qualifiers.Read("qualifiers.dat");
	for(int iqual=0;iqual<qualifiers.nqualifiers;iqual++){
		dnchdy=0;
		ncolls=0;
		npartstot=0;
		b3d->SetQualifier(qualifiers.qualifier[iqual]);
		qualifiers.SetPars(&(b3d->parmap),iqual);
		nevents=parameter::getI(b3d->parmap,"B3D_NEVENTSMAX",10);
		printf("nevents=%d\n",nevents);
		b3d->sampler->ReadVolumeElements2D();
		printf("check\n");
		for(ievent=1;ievent<=nevents;ievent++){
			b3d->Reset();
			b3d->randy->reset(-ievent);
			nparts=b3d->sampler->MakeB3DEvent();
			b3d->PerformAllActions();
			dnchdy+=b3d->WriteOSCAR(ievent);
			printf("nscatter=%g, nmerges=%g, ndecays=%g,  ncellexits=%g, nregenerate=%g\n",
			double(b3d->nscatter),double(b3d->nmerge),double(b3d->ndecay),double(b3d->nexit),double(b3d->nregenerate));
			ncolls+=b3d->nscatter+b3d->nmerge;
			nbaryons+=b3d->nbaryons;
			nannihilate+=b3d->nannihilate;
			nregen+=b3d->nregenerate;
			npartstot+=b3d->FinalPartMap.size();
      printf("###### finished event %d ####### %d FS parts #########\n",ievent,int(b3d->FinalPartMap.size()));
		}
		norm=nevents*b3d->NSAMPLE;
		printf("<# of FS parts>=%g, <# nB>=%g, <# collisions>=%g, <# annihilations>=%g <# regenerations=%g>\n",double(npartstot)/norm,double(nbaryons)/norm,double(ncolls)/norm,double(nannihilate)/norm,double(nregen)/norm);
	}
	return 0;
}
