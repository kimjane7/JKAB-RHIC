#ifndef __RESONANCES_H__
#define __RESONANCES_H__

#include "b3d.h"

class CResInfo;
class CBranchInfo;

typedef unordered_map<long int,CResInfo *> CResInfoMap;
typedef pair<long int, CResInfo*> CResInfoPair;
typedef vector<CBranchInfo *> CBranchList; //gives branchlist name

class CB3D;
class CBranchInfo{
public:
	vector<CResInfo *> resinfoptr; //pointers for resinfo
	double branching; 
	CBranchInfo();
};

class CMerge{
public:
	CMerge(CResInfo *resinfo,double branching, int L);
	CResInfo *resinfo;
	int L;
	double branching;
	CMerge *next;
};

class CResInfo{
public:
	int ires;
	double mass;
	double spin;
	double width; 
	double minmass;
	string name;
	int code;
	int charge;
	int strange;
	int baryon;	
	int G_Parity;
	bool decay; //false if stable, true if can decay. check if true
	CBranchList branchlist; 
	CBranchInfo	*bptr_minmass;
	void Print();
	void DecayGetResInfoPtr(int &nbodies,array<CResInfo *,5> &daughterresinfo);
	void DecayGetResInfoPtr_minmass(int &nbodies,array<CResInfo *,5> &daughterresinfo);
	bool CheckForDaughters(int code);
	bool CheckForNeutral();
	double GenerateThermalMass(double maxweight, double T);
	CResInfo();
	static CRandom *ranptr;
};

class CResList{
public:
	CResList();
	~CResList();
	CResList(parameterMap* parmap_in);
	CResInfoMap resmap;
	CResInfo *GetResInfoPtr(int ID);
	void ReadResInfo();
	void CalcEoS(double T,double &epsilon,double &P,double &nhadrons,vector<double> &density,vector<double> &boseweight,vector<double> &maxweight);
	void CalcEoS(double T,double &epsilon,double &P,double &nhadrons,double &cs2,vector<double> &density);
	void CalcEoS(double T_0,double T_F,double delT);
	void CalcEoSandChi(double T);
	void CalcEosandKubo(double T,double &epsilon,double &P,double &nhadrons,double &sdens,double &kubo);
	void freegascalc_onespecies(double m,double T,double &e,double &p,double &dens,double &sigma2,double &dedt);
	void freegascalc_onespecies_finitewidth(double m, double m1, double m2, double T,double width,double minmass,double &epsilon,double &P,double &dens,double &sigma2,double &dedt, double &maxweight);
	parameterMap *parmap;
	CMerge ***MergeArray;
	double **SigmaMaxArray;
	static CB3D *b3d;
	//void freegascalc_onespecies_offshell(CResInfo *resinfo,double T,double &epsilon,double &P,double &dens,double &sigma2,double &dedt);
};




#endif
