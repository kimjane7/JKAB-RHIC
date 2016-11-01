#include "CHydro.h"
#include <fstream>
#include "hydroDef.h"
#include <cstdlib>
#include <cstdio> 
#include <cmath>

#include <ctime>
#include <math.h>
#include <string.h>
#include "CMesh.h" 
#include <fullMesh.h>
#include <cornelius.h>
	
void CHydro::FindTeddySurface(){
	CMesh *newMesh=tempMesh;
  CMesh *prevMesh=offMesh;
	double newtau=newMesh->getTau();
	double prevtau=prevMesh->getTau();

	double omega[4];
	bool allabove=true,allbelow=true;
	double DX=newMesh->mDx,DETA=newMesh->mDn;
	double delx[4], GradTdotomega;
	double lowX1=0, lowY1=0, lowN1=0;
	int ix,iy,iz,itau,idx,idy,idz,idtau,i,x,y,z;
	double T[2][2][2][2];
	double Tplus[4],Tminus[4],GradT[4];
	double taubar=0.5*(newtau+prevtau);
	delx[0]=newtau-prevtau;
	delx[1]=delx[2]=newMesh->mDn;
	delx[3]=newMesh->mDn*taubar;
	double Omega[4], GradT2,dV;
	double Xbar,Ybar,Nbar, E, P, Temp, Pixy;
	
	for(int iz=lowN1; iz<newMesh->getNSize(); iz++){
		for(int iy=lowY1; iy<newMesh->getYSize(); iy++){
			for(int ix=lowX1; ix<newMesh->getXSize(); ix++){

				newMesh->getXSize(),newMesh->getYSize(),newMesh->getNSize(),
				prevMesh->getXSize(),prevMesh->getYSize(),prevMesh->getNSize();
				for(idtau=0;idtau<2;idtau++){
					for(idx=0;idx<2;idx++){
						for(idy=0;idy<2;idy++){
							for(idz=0;idz<2;idz++){
								if(idtau==0)
									T[idtau][idx][idy][idz]=prevMesh->getT(iz+idz,ix+idx,iy+idy);
								if(idtau==1)
									T[idtau][idx][idy][idz]=newMesh->getT(iz+idz,ix+idx,iy+idy);
								if(T[idtau][idx][idy][idz]>mFoTemp)
									allbelow=false;
								if(T[idtau][idx][idy][idz]<mFoTemp)
									allabove=false;
							}
						}
					}
				}
				for(idtau=0;idtau<2;idtau++){
					for(idx=0;idx<2;idx++){
						for(idy=0;idy<2;idy++){
							for(idz=0;idz<2;idz++){
								Xbar = ((newMesh->getCell(iz,ix,iy))->getX() + (newMesh->getCell(iz,ix+idx,iy))->getX())*0.5;
							  Ybar = ((newMesh->getCell(iz,ix,iy))->getY() + (newMesh->getCell(iz,ix,iy+idy))->getY())*0.5;
							  Nbar = ((newMesh->getCell(iz,ix,iy))->getZ() + (newMesh->getCell(iz+idz,ix,iy))->getZ())*0.5;
								Temp = ((newMesh->getCell(iz,ix,iy))->getT() + (newMesh->getCell(iz+idz,ix+idx,iy+idy))->getZ())*0.5;
								E = ((newMesh->getCell(iz,ix,iy))->getE()+ (newMesh->getCell(iz+idz,ix+idx,iy+idy))->getE())*0.5;							
								P = ((newMesh->getCell(iz,ix,iy))->getP() + (newMesh->getCell(iz+idz,ix+idx,iy+idy))->getP())*0.5;
								Pixy = ((newMesh->getCell(iz,ix,iy))->getPixy(ix,iy) + (newMesh->getCell(iz+idz,ix+idx,iy+idy))->getPixy(ix,iy))*0.5;
								
							}
						}
					}
				
				if(allabove==false && allbelow==false){
					Tplus[0]=Tminus[0]=Tplus[1]=Tminus[1]=Tplus[2]=Tminus[2]=Tplus[3]=Tminus[3];
					for(idtau=0;idtau<2;idtau++){
						for(idx=0;idx<2;idx++){
							for(idy=0;idy<2;idy++){
								for(idz=0;idz<2;idz++){									
									if(idtau==0)
										Tminus[0]+=0.125*T[idtau][idx][idy][idz];
									else
										Tplus[0]+=0.125*T[idtau][idx][idy][idz];
									if(idx==0)
										Tminus[1]+=0.125*T[idtau][idx][idy][idz];
									else
										Tplus[1]+=0.125*T[idtau][idx][idy][idz];
									if(idy==0)
										Tminus[2]+=0.125*T[idtau][idx][idy][idz];
									else
										Tplus[2]+=0.125*T[idtau][idx][idy][idz];
									if(idz==0)
										Tminus[3]+=0.125*T[idtau][idx][idy][idz];
									else
										Tplus[3]+=0.125*T[idtau][idx][idy][idz];
								}
							}
						}
					}
					allabove=allbelow=true;
					for(i=0;i<4;i++){
						if(Tplus[i]>mFoTemp)
							allbelow=false;
						else
							allabove=false;
					}
					if(allabove==false && allbelow==false){
						for(i=0;i<4;i++)
							GradT[i]=(Tplus[i]-Tminus[i])/delx[i];
						omega[0]=omega[1]=omega[2]=omega[3]=0.0;
						for(i=0;i<4;i++){
							if((Tplus[i]>mFoTemp && Tminus[i]<mFoTemp) || (Tplus[i]<mFoTemp && Tminus[i]>mFoTemp)){
								omega[i]+=(dV/delx[i])*(Tplus[i]-Tminus[i])*fabs(Tplus[i]-Tminus[i]);
							}
						}
						GradTdotomega=omega[0]*GradT[0]+omega[1]*GradT[1]+omega[2]*GradT[2]+omega[3]*GradT[3];
						GradT2=GradT[0]*GradT[0]+GradT[1]*GradT[1]+GradT[2]*GradT[2]+GradT[3]*GradT[3];
						for(i=0;i<4;i++){
							Omega[i]=GradT[i]*(GradTdotomega/GradT2); 
						  
	            /*&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
	              &                                                                                                                            &
		            &   Output order: tau, x, y, z, ux, uy, uz, omega0, omegax, omegay, omegaz, pixx, pixy, pixz, piyy, piyz, pizz, T, E, P      &
		            &					                                                                                                                   &
	              &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*/
	
	
							fprintf(fTeddyHyper,"%f, %f, %f, %f, %g, %g, %g, %g, %g, %g, %g, %g  \n", taubar, Xbar, Ybar, Nbar, Omega[0],Omega[1],Omega[2],Omega[3],Pixy,Temp,E,P);		
						  printf(" \n ");
						}
					}
				}
			}
		}
	}
}
}