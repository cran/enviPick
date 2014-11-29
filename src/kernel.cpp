/*
Partitioning, clustering & peak detection for LC-MS centroided data
author: Martin Loos, Martin.Loos@eawag.ch
Copyright (c) 2014 Eawag. All rights reserved.
*/

#include <list>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <vector>
#include "auxil.h"

using namespace std;

extern "C"{

      /************************************************************************/
      /* agglomerative partitioning *******************************************/
      /************************************************************************/
      SEXP agglom(     SEXP mz, /* must be sorted */
                       SEXP rt,
                       SEXP ppm,
                       SEXP dmz,
                       SEXP drt
                      ){


            PROTECT(mz = AS_NUMERIC(mz));
            PROTECT(rt = AS_NUMERIC(rt));
            PROTECT(ppm = AS_INTEGER(ppm));
            PROTECT(dmz = AS_NUMERIC(dmz));
            PROTECT(drt = AS_NUMERIC(drt));
            double *mass;
            mass = NUMERIC_POINTER(mz);
            double *ret;
            ret = NUMERIC_POINTER(rt);
            int ppm2 = INTEGER_VALUE(ppm);
            double dmass = NUMERIC_VALUE(dmz);
            double dret = NUMERIC_VALUE(drt);
            int n,m,p,k=0,found=0;
            double lowmass, highmass, lowret, highret;
            int leng = LENGTH(rt);
            SEXP outit;
            PROTECT(outit = NEW_INTEGER(leng));
            int *at;
            at = INTEGER_POINTER(outit);
            for(n=0;n<leng;n++){*(at+n) = 0;}
            SETLENGTH(outit,leng);
            int *these;
            these = new int[leng];
            int *those;
            those = new int[leng];
            int untilthese=0,untilthose=0,atthese;

            for(n=0;n<leng;n++){
                if(*(at+n) == 0){ /* unassigned entry ? */
                    k++;
                    *(at+n)=k;
                    found=1;
                    these[0]=n;
                    untilthese=1;
                    atthese=1;
                    while(found == 1){
                        found=0;
                        if(atthese==1){
                            untilthose=0;
                            for(m=0;m<untilthese;m++){
                                if(ppm2 == 1){
                                   lowmass=(*(mass+these[m])-((*(mass+these[m])*dmass)/1E6));
                                   highmass=(*(mass+these[m])+((*(mass+these[m])*dmass)/1E6));
                                }else{
                                   lowmass=(*(mass+these[m])-(dmass));
                                   highmass=(*(mass+these[m])+(dmass));
                                }
                                lowret=(*(ret+these[m])-dret);
                                highret=(*(ret+these[m])+dret);
                                if(these[m]>0){
                                    for(p=(these[m]-1);p>=0;p--){ /* towards lower mass */
                                        if((*(mass+p))<=(lowmass)){
                                            break;
                                        }else{
                                            if(*(at+p)==0){
                                                if((*(ret+p)>=lowret)&&(*(ret+p)<=highret)){
                                                    those[untilthose]=p;
                                                    untilthose++;
                                                    *(at+p)=k;
                                                    found=1;
                                                }
                                            }
                                        }
                                    }
                                }
                                if(these[m]<(leng-1)){
                                    for(p=(these[m]+1);p<leng;p++){ /* towards higher mass */
                                        if((*(mass+p))>=(highmass)){
                                            break;
                                        }else{
                                            if(*(at+p)==0){
                                                if((*(ret+p)>=lowret)&&(*(ret+p)<=highret)){
                                                    those[untilthose]=p;
                                                    untilthose++;
                                                    *(at+p)=k;
                                                    found=1;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                            atthese=0;
                        }else{
                            untilthese=0;
                            for(m=0;m<untilthose;m++){
                                if(ppm2 == 1){
                                   lowmass=(*(mass+those[m])-((*(mass+those[m])*dmass)/1E6));
                                   highmass=(*(mass+those[m])+((*(mass+those[m])*dmass)/1E6));
                                }else{
                                   lowmass=(*(mass+those[m])-(dmass));
                                   highmass=(*(mass+those[m])+(dmass));
                                }
                                lowret=(*(ret+those[m])-dret);
                                highret=(*(ret+those[m])+dret);
                                if(those[m]>0){
                                    for(p=(those[m]-1);p>=0;p--){
                                        if((*(mass+p))<=(lowmass)){
                                            break;
                                        }else{
                                            if(*(at+p)==0){
                                                if((*(ret+p)>=lowret)&(*(ret+p)<=highret)){
                                                    these[untilthese]=p;
                                                    untilthese++;
                                                    *(at+p)=k;
                                                    found=1;
                                                }
                                            }
                                        }
                                    }
                                }
                                if(those[m]<(leng-1)){
                                    for(p=(those[m]+1);p<leng;p++){
                                       if((*(mass+p))>=(highmass)){
                                           break;
                                        }else{
                                             if(*(at+p)==0){
                                                 if((*(ret+p)>=lowret)&(*(ret+p)<=highret)){
                                                     these[untilthese]=p;
                                                     untilthese++;
                                                     *(at+p)=k;
                                                     found=1;
                                                 }
                                             }
                                        }
                                    }
                                }
                            }
                            atthese=1;
                        }
                    }
                }
            }

            delete[] these;
            delete[] those;
            SETLENGTH(outit,leng);
            UNPROTECT(6);
            return outit;

      }

      /************************************************************************/
      /* assemble indices for return value of agglom **************************/
      /************************************************************************/
      SEXP indexed(      SEXP index, /* must be sorted */
				         SEXP intensity,
                         SEXP minpeak,
                         SEXP maxint,
                         SEXP maxindex
                         ){

            PROTECT(index = AS_NUMERIC(index));
            PROTECT(intensity = AS_NUMERIC(intensity));
            PROTECT(minpeak = AS_INTEGER(minpeak));
            PROTECT(maxint = AS_NUMERIC(maxint));
            PROTECT(maxindex = AS_NUMERIC(maxindex));
            double *ind;
            ind = NUMERIC_POINTER(index);
            double *inte;
            inte = NUMERIC_POINTER(intensity);
            int minpeak2 = INTEGER_VALUE(minpeak);
            double maxint2 = NUMERIC_VALUE(maxint);
			double tempmax=0;
            int maxind = INTEGER_VALUE(maxindex);
			int leng = LENGTH(index);
            int n,from,to,counted,atind;
            SEXP outit;
            PROTECT(outit = allocMatrix(INTSXP, maxind, 3));
            int *at;
            at = INTEGER_POINTER(outit);
            for(n=0;n<(maxind*3);n++){
                *(at+n) = 0;
            }

            from=1;
            to=1;
            counted=1;
            tempmax=*inte;
            atind=0;
            for(n=1;n<leng;n++){
                if(*(ind+n)!=*(ind+n-1)){
                	if( ((tempmax>=maxint2) || (counted>=minpeak2)) && (*(ind+n-1)!=0) ){
                   		*(at+atind)=from;
                        *(at+maxind+atind)=to;
                   		*(at+(maxind*2)+atind)=counted;
                   		atind++;
                    };
                	from=(n+1);
                	to=from;
                    counted=1;
                    tempmax=*(inte+n);
                }else{
                	if(*(inte+n)>tempmax){
               	        tempmax=*(inte+n);
                	};
                    to++;
                    counted++;
                }
            }
            n--;
            if(  ((tempmax>=maxint2)||(counted>=minpeak2)) && (*(ind+n)!=0)  ){
                *(at+atind)=from;
                *(at+maxind+atind)=to;
                *(at+(maxind*2)+atind)=counted;
            }

            SETLENGTH(outit,maxind*3);
            UNPROTECT(6);
            return outit;

      }


      /************************************************************************/
      /* assemble agglom part ID vector for measurements **********************/
      /************************************************************************/
      SEXP partID(       SEXP index,
                         SEXP leng
                         ){

            PROTECT(index = AS_INTEGER(index));
            PROTECT(leng = AS_INTEGER(leng));
            int lengit = LENGTH(index);
            lengit=lengit/3;
            int *indexed;
            indexed = INTEGER_POINTER(index);
            int leng2 = INTEGER_VALUE(leng);
            SEXP outit;
            PROTECT(outit = NEW_INTEGER(leng2));
            int *at;
            int n,m,id=1;
            at = INTEGER_POINTER(outit);
            for(n=0;n<leng2;n++){
                *(at+n) = 0;
            }

            for(n=0;n<lengit;n++){
                for((m=*(indexed+n));(m<=*(indexed+lengit+n));m++){
                    *(at+(m-1))=id;
                }
                id++;
            }

            UNPROTECT(3);
            return outit;
      }


      /************************************************************************/
      /* m/z partitioning *****************************************************/
      /************************************************************************/
       SEXP masspart(  SEXP mz, /* must be sorted */
                       SEXP dmz,
                       SEXP ppm
                      ){

            PROTECT(mz = AS_NUMERIC(mz));
            PROTECT(dmz = AS_NUMERIC(dmz));
            PROTECT(ppm = AS_INTEGER(ppm));

            double diffmass;
            double *mz2;
            mz2 = NUMERIC_POINTER(mz);
            double dmz2 = NUMERIC_VALUE(dmz);
            double ppm2 = INTEGER_VALUE(ppm);
            int n,m;
            int leng = LENGTH(mz);
            int *peaked;
            peaked = new int[leng];

            peaked[0]=1;
            m=1;
            if(leng>1){
                for(n=1;n<leng;n++){
                    if(ppm2 == 1){
                        diffmass=((fabs(*(mz2+n)-*(mz2+n-1)))/ *(mz2+n) *1E6);
                    }else{
                        diffmass=fabs(*(mz2+n)-*(mz2+n-1));
                    }
                    if(diffmass>dmz2){
                        m++;
                        peaked[n]=m;
                    }else{
                        peaked[n]=m;
                    }
                }
            }

            SEXP outit;
            PROTECT(outit = NEW_INTEGER(leng));
            int *at;
            at = INTEGER_POINTER(outit);
            for(n=0;n<leng;n++){
                   *(at+n) = peaked[n];
            }
            SETLENGTH(outit,leng);
            UNPROTECT(4);
            delete[] peaked;
            return outit;

       }


      /************************************************************************/
      /* RT partitioning ******************************************************/
      /************************************************************************/
        SEXP rtpart(  SEXP rt, /*must be sorted*/
                      SEXP drt
                      ){

            PROTECT(rt = AS_NUMERIC(rt));
            PROTECT(drt = AS_NUMERIC(drt));

            double *rt2;
            rt2 = NUMERIC_POINTER(rt);
            double drt2 = NUMERIC_VALUE(drt);
            int n,m;
            int leng = LENGTH(rt);
            int *peaked;
            peaked = new int[leng];

            peaked[0]=1;
            m=1;
            if(leng>1){
                for(n=1;n<leng;n++){
                    if(fabs(*(rt2+n)-*(rt2+n-1))>drt2){
                        m++;
                        peaked[n]=m;
                    }else{
                        peaked[n]=m;
                    }
                }
            }

            SEXP outit;
            PROTECT(outit = NEW_INTEGER(leng));
            int *at;
            at = INTEGER_POINTER(outit);
            for(n=0;n<leng;n++){
                   *(at+n) = peaked[n];
            }
            SETLENGTH(outit,leng);
            UNPROTECT(3);
            delete[] peaked;
            return outit;

       }


      /************************************************************************/
      /* EIC clustering *******************************************************/
      /************************************************************************/
      SEXP getEIC(      SEXP mz,
                        SEXP RT,
                        SEXP intens,
                        SEXP orderedint,
                        SEXP orderedret,
                        SEXP dmzdens,
                        SEXP ppm2,
                        SEXP drtdens,
                        SEXP merged2
                        ){

           PROTECT(mz = AS_NUMERIC(mz));
           PROTECT(RT = AS_NUMERIC(RT));
           PROTECT(intens = AS_NUMERIC(intens));
           PROTECT(orderedint = AS_INTEGER(orderedint));
           PROTECT(orderedret = AS_INTEGER(orderedret));
           PROTECT(dmzdens = AS_NUMERIC(dmzdens));
           PROTECT(ppm2 = AS_INTEGER(ppm2));
           PROTECT(drtdens = AS_NUMERIC(drtdens));
           PROTECT(merged2 = AS_INTEGER(merged2));
           double *ret, *mass, *intensity;
           mass = NUMERIC_POINTER(mz);
           ret = NUMERIC_POINTER(RT);
           intensity = NUMERIC_POINTER(intens);
           int *ordint, *ordret;
           ordint = INTEGER_POINTER(orderedint);
           ordret = INTEGER_POINTER(orderedret);
           double dmzdens2 = NUMERIC_VALUE(dmzdens);
           int ppm3 = INTEGER_VALUE(ppm2);
           double drtdens2 = NUMERIC_VALUE(drtdens);
           int merged3 = INTEGER_VALUE(merged2);
           int leng = LENGTH(RT);
           int m,n,i,k,clustnumb,maxat=0,maxit=0;
           double delmz;
           SEXP clusters;
           PROTECT(clusters = allocMatrix(REALSXP, leng, 13));
           double *clus;
           clus = REAL(clusters);
           for(m=0;m<13;m++){
               for(n=0;n<leng;n++){
                   clus[(m*leng)+n]=0;
               }
           }
           int *at;
           at = new int[leng];

           /* initialize with most intense measurement ************************/
           clustnumb=1;
           if(ppm3==1){delmz=((dmzdens2**(mass+(*(ordint)-1)))/1e6);}else{delmz=dmzdens2;}
           clus[0]=(*(mass+(*(ordint)-1))-(2*delmz));       /* low mass boundary **************/
           clus[(1*leng)]=(*(mass+(*(ordint)-1))+(2*delmz));/* high mass boundary *************/
           clus[(2*leng)]=(*(ret+(*(ordint)-1))-drtdens2);  /* low RT boundary ****************/
           clus[(3*leng)]=(*(ret+(*(ordint)-1))+drtdens2);  /* high RT boundary ***************/
           clus[(4*leng)]=1;                                /* number of measurements *********/
           clus[(5*leng)]=*(mass+(*(ordint)-1));            /* mass sum ***********************/
           clus[(6*leng)]=clustnumb;                        /* cluster ID *********************/
           clus[(7*leng)]=0;                                /* merged (1) or not (0)? *********/
           clus[(8*leng)]=*(intensity+(*(ordint)-1));       /* maximum intensity in a cluster */
           clus[(9*leng)+(*(ordint)-1)]=clustnumb;          /* cluster ID for measurement *****/
           clus[(10*leng)]=0;                               /* variance, set later if merged **/
           clus[(11*leng)]=*(mass+(*(ordint)-1));           /* lowest mass in cluster *********/
           clus[(12*leng)]=*(mass+(*(ordint)-1));           /* highest mass in cluster ********/

           /* assign all other peaks ******************************************/
           for(n=1;n<leng;n++){
               /* check for possible fit to existing clusters *****************/
               maxat=0;
               for(m=0;m<clustnumb;m++){
                if(*(mass+(*(ordint+n)-1))<=clus[(1*leng)+m]){
                   if(*(mass+(*(ordint+n)-1))>=clus[(0*leng)+m]){
                           if(*(ret+(*(ordint+n)-1))>=clus[(2*leng)+m]){
                               if(*(ret+(*(ordint+n)-1))<=clus[(3*leng)+m]){
                                       at[maxat]=m;
                                       maxat++;
                               }
                           }
                       }
                   }
               }
               /* (a) not assignable - create new cluster *********************/
               if(maxat==0){
                   clustnumb++;
                   clus[(9*leng)+(*(ordint+n)-1)]=clustnumb;
                   if(ppm3==1){delmz=((dmzdens2**(mass+(*(ordint+n)-1)))/1e6);}else{delmz=dmzdens2;}
                   clus[0+(clustnumb-1)]=(*(mass+(*(ordint+n)-1))-(2*delmz));
                   clus[(1*leng)+(clustnumb-1)]=(*(mass+(*(ordint+n)-1))+(2*delmz));
                   clus[(2*leng)+(clustnumb-1)]=(*(ret+(*(ordint+n)-1))-drtdens2);
                   clus[(3*leng)+(clustnumb-1)]=(*(ret+(*(ordint+n)-1))+drtdens2);
                   clus[(4*leng)+(clustnumb-1)]=1;
                   clus[(5*leng)+(clustnumb-1)]=*(mass+(*(ordint+n)-1));
                   clus[(6*leng)+(clustnumb-1)]=clustnumb;
                   clus[(7*leng)+(clustnumb-1)]=0;
                   clus[(8*leng)+(clustnumb-1)]=*(intensity+(*(ordint+n)-1));
                   clus[(11*leng)+(clustnumb-1)]=*(mass+(*(ordint+n)-1));
                   clus[(12*leng)+(clustnumb-1)]=*(mass+(*(ordint+n)-1));
                   continue;
               }
               /* check for other peaks at same RT ****************************/
               maxit=maxat;
               for(i=0;i<leng;i++){ /* find RT-order index */
                   if(*(ordint+n)==*(ordret+i)){
                       k=i;
                       break;
                   }
               }

               if(k>0){ /* backward over RT-order */
                   for(i=(k-1);i>=0;i--){
                       if(  *(ret+(*(ordret+k)-1))==*(ret+(*(ordret+i)-1))  ){
                           if((clus[(9*leng)+(*(ordret+i)-1)])==0){
                               continue;
                           }else{
                               for(m=0;m<maxat;m++){
                                   if(clus[(9*leng)+(*(ordret+i)-1)]==(at[m]+1)){
                                       at[m]=-9999;
                                       maxit--;
                                       break;
                                   }
                               }
                           }
                       }else{
                           break;
                       }
                   }
               }

               if(k<(leng-1)){ /* forward over RT-order */
                   for(i=(k+1);i<leng;i++){
                       if(*(ret+(*(ordret+k)-1))==*(ret+(*(ordret+i)-1))){
                           if((clus[(9*leng)+(*(ordret+i)-1)])==0){
                               continue;
                           }else{
                               for(m=0;m<maxat;m++){
                                   if(clus[(9*leng)+(*(ordret+i)-1)]==(at[m]+1)){
                                       at[m]=-9999;
                                       maxit--;
                                       break;
                                   }
                               }
                           }
                       }else{
                           break;
                       }
                   }
               }

               if(maxit<maxat){
                   i=maxat;
                   k=0;
                   for(m=0;m<i;m++){
                       if(at[m]==-9999){
                           k++;
                           maxat--;
                       }else{
                           at[m-k]=at[m];
                       }
                   }
               }
               maxat=maxit;
               /* (a) not assignable - create new cluster *********************/
               if(maxat==0){
                   clustnumb++;
                   clus[(9*leng)+(*(ordint+n)-1)]=clustnumb;
                   if(ppm3==1){delmz=((dmzdens2**(mass+(*(ordint+n)-1)))/1e6);}else{delmz=dmzdens2;}
                   clus[0+(clustnumb-1)]=(*(mass+(*(ordint+n)-1))-(2*delmz));
                   clus[(1*leng)+(clustnumb-1)]=(*(mass+(*(ordint+n)-1))+(2*delmz));
                   clus[(2*leng)+(clustnumb-1)]=(*(ret+(*(ordint+n)-1))-drtdens2);
                   clus[(3*leng)+(clustnumb-1)]=(*(ret+(*(ordint+n)-1))+drtdens2);
                   clus[(4*leng)+(clustnumb-1)]=1;
                   clus[(5*leng)+(clustnumb-1)]=*(mass+(*(ordint+n)-1));
                   clus[(6*leng)+(clustnumb-1)]=clustnumb;
                   clus[(7*leng)+(clustnumb-1)]=0;
                   clus[(8*leng)+(clustnumb-1)]=*(intensity+(*(ordint+n)-1));
                   clus[(11*leng)+(clustnumb-1)]=*(mass+(*(ordint+n)-1));
                   clus[(12*leng)+(clustnumb-1)]=*(mass+(*(ordint+n)-1));
                   continue;
               }
               /* (b) fits exactly to one cluster *****************************/
               if(maxat==1){
                   clus[(9*leng)+(*(ordint+n)-1)]=(at[0]+1);
                   if(ppm3==1){delmz=((dmzdens2**(mass+(*(ordint+n)-1)))/1e6);}else{delmz=dmzdens2;}
                   /* shrink lower & upper mass bounds */
                   if(clus[(0*leng)+(at[0])]<(*(mass+(*(ordint+n)-1))-(2*delmz))){clus[(0*leng)+(at[0])]=(*(mass+(*(ordint+n)-1))-(2*delmz));}
                   if(clus[(1*leng)+(at[0])]>(*(mass+(*(ordint+n)-1))+(2*delmz))){clus[(1*leng)+(at[0])]=(*(mass+(*(ordint+n)-1))+(2*delmz));}
                   /* adapt lower & upper RT bounds */
                   if(clus[(2*leng)+(at[0])]>(*(ret+(*(ordint+n)-1))-drtdens2)){clus[(2*leng)+(at[0])]=(*(ret+(*(ordint+n)-1))-drtdens2);}
                   if(clus[(3*leng)+(at[0])]<(*(ret+(*(ordint+n)-1))+drtdens2)){clus[(3*leng)+(at[0])]=(*(ret+(*(ordint+n)-1))+drtdens2);}
                   clus[(4*leng)+(at[0])]=clus[(4*leng)+(at[0])]+1;
                   clus[(5*leng)+(at[0])]=clus[(5*leng)+(at[0])]+*(mass+(*(ordint+n)-1));
                   if(*(mass+(*(ordint+n)-1))<clus[(11*leng)+(at[0])]){clus[(11*leng)+(at[0])]=*(mass+(*(ordint+n)-1));};
                   if(*(mass+(*(ordint+n)-1))>clus[(12*leng)+(at[0])]){clus[(12*leng)+(at[0])]=*(mass+(*(ordint+n)-1));};
                   continue;
               }
               /* (c) fits to several clusters - m/z density clustering *******/
               if(maxat>1){
                   delmz=fabs(*(mass+(*(ordint+n)-1))- ( clus[(5*leng)+(at[0])]/clus[(4*leng)+(at[0])] ));
                   for(m=0;m<maxat;m++){
                       if(
                           fabs(*(mass+(*(ordint+n)-1))- ( clus[(5*leng)+(at[m])]/clus[(4*leng)+(at[m])] ))<
                           delmz
                       ){
                           at[0]=at[m];
                           delmz=fabs(*(mass+(*(ordint+n)-1))- ( clus[(5*leng)+(at[0])]/clus[(4*leng)+(at[0])] ));
                       }
                   }
                   clus[(9*leng)+(*(ordint+n)-1)]=(at[0]+1);
                   if(ppm3==1){delmz=((dmzdens2**(mass+(*(ordint+n)-1)))/1e6);}else{delmz=dmzdens2;}
                   /* shrink lower & upper mass bounds */
                   if(clus[(0*leng)+(at[0])]<(*(mass+(*(ordint+n)-1))-(2*delmz))){clus[(0*leng)+(at[0])]=(*(mass+(*(ordint+n)-1))-(2*delmz));}
                   if(clus[(1*leng)+(at[0])]>(*(mass+(*(ordint+n)-1))+(2*delmz))){clus[(1*leng)+(at[0])]=(*(mass+(*(ordint+n)-1))+(2*delmz));}
                   /* adapt lower & upper RT bounds */
                   if(clus[(2*leng)+(at[0])]>(*(ret+(*(ordint+n)-1))-drtdens2)){clus[(2*leng)+(at[0])]=(*(ret+(*(ordint+n)-1))-drtdens2);}
                   if(clus[(3*leng)+(at[0])]<(*(ret+(*(ordint+n)-1))+drtdens2)){clus[(3*leng)+(at[0])]=(*(ret+(*(ordint+n)-1))+drtdens2);}
                   clus[(4*leng)+(at[0])]=clus[(4*leng)+(at[0])]+1;
                   clus[(5*leng)+(at[0])]=clus[(5*leng)+(at[0])]+*(mass+(*(ordint+n)-1));
                   if(*(mass+(*(ordint+n)-1))<clus[(11*leng)+(at[0])]){clus[(11*leng)+(at[0])]=*(mass+(*(ordint+n)-1));};
                   if(*(mass+(*(ordint+n)-1))>clus[(12*leng)+(at[0])]){clus[(12*leng)+(at[0])]=*(mass+(*(ordint+n)-1));};
                   continue;
              }
           }

           /* *************************************************************** */
           /* merge isobaric mass cluster *************************************/
           if((merged3==1)&(clustnumb>1)){

               /* set mean ****************************************************/
               for(n=0;n<clustnumb;n++){
                   clus[(5*leng)+(n)]=(clus[(5*leng)+(n)]/clus[(4*leng)+(n)]);
               }
               /* define exclusion-list of cluster with same RT ***************/
               std::vector<int> blackone;
               std::vector<int> blacktwo;
               int blacksize;
               for(n=1;n<leng;n++){
                   for(k=(n-1);k>=0;k--){
                       if(*(ret+(*(ordret+k)-1))==*(ret+(*(ordret+n)-1))){
                           if( /* confine to overlapping masses */
                                (clus[(0*leng)+(int(clus[(9*leng)+(*(ordret+n)-1)])-1)] < clus[(1*leng)+(int(clus[(9*leng)+(*(ordret+k)-1)])-1)]) &
                                (clus[(1*leng)+(int(clus[(9*leng)+(*(ordret+n)-1)])-1)] > clus[(0*leng)+(int(clus[(9*leng)+(*(ordret+k)-1)])-1)])

                           ){          /* smallest index always first */
                               if((clus[(9*leng)+(*(ordret+k)-1)])<(clus[(9*leng)+(*(ordret+n)-1)])){
                                       blackone.push_back(int(clus[(9*leng)+(*(ordret+k)-1)]));
                                       blacktwo.push_back(int(clus[(9*leng)+(*(ordret+n)-1)]));
                               }else{
                                       blackone.push_back(int(clus[(9*leng)+(*(ordret+n)-1)]));
                                       blacktwo.push_back(int(clus[(9*leng)+(*(ordret+k)-1)]));
                               }
                           }
                       }else{
                           break;
                       }
                   }
               }
               blacksize=blackone.size();
               /* find mergeable clusters: indices & mass differences *********/
               /* cross-check for overlaps in RT between cluster **************/
               std::vector<int> clusterone;
               std::vector<int> clustertwo;
               std::vector<double> clusterdiff;
               int doit;
               for(n=0;n<(clustnumb-1);n++){
                   for(m=(n+1);m<clustnumb;m++){
                       if(
                           (clus[(0*leng)+(n)]<clus[(1*leng)+(m)]) &
                           (clus[(1*leng)+(n)]>clus[(0*leng)+(m)])
                       ){ /* check for overlapping masses */
                           if(
                               (
                                   (clus[(0*leng)+(n)]<=clus[(11*leng)+(m)]) &
                                   (clus[(1*leng)+(n)]>=clus[(12*leng)+(m)])
                                )||(
                                   (clus[(0*leng)+(m)]<=clus[(11*leng)+(n)]) &
                                   (clus[(1*leng)+(m)]>=clus[(12*leng)+(n)])
                                )
                           ){
                               if(blacksize==0){
                                   clusterone.push_back(n+1);
                                   clustertwo.push_back(m+1);
                                   clusterdiff.push_back(fabs(clus[(5*leng)+(n)]-clus[(5*leng)+(m)]));
                               }else{
                                   doit=0;
                                   for(k=0;k<blacksize;k++){
                                       if(blackone[k]==(n+1)){
                                            if(blacktwo[k]==(m+1)){
                                                doit=1;
                                                break;
                                            }
                                       }
                                   }
                                   if(doit==0){
                                       clusterone.push_back(n+1);
                                       clustertwo.push_back(m+1);
                                       clusterdiff.push_back(fabs(clus[(5*leng)+(n)]-clus[(5*leng)+(m)]));
                                    }
                               }
                           }
                       }
                   }
                }
                /* merge cluster ***********************************************/
                int mergesize=clusterone.size(), stay, gone;
                std::vector<int> clusterase;
               if(mergesize>0){
                   while(mergesize>0){
                        /* find cluster pair closest in mean mass difference */
                        delmz=clusterdiff[0];
                        m=0;
						if(mergesize>1){
                            for(k=1;k<mergesize;k++){
                               if(clusterdiff[k]<delmz){
                                   m=k;
                                   delmz=clusterdiff[k];
                               }
                            }
                    	}
                        /* reassign measurements */
						for(n=0;n<leng;n++){
                            if(clus[(9*leng)+n]==clustertwo[m]){
                                clus[(9*leng)+n]=clusterone[m];
                            }
                        }
						/* update cluster entries */
                        clus[(4*leng)+(clusterone[m]-1)]=(clus[(4*leng)+(clusterone[m]-1)]+clus[(4*leng)+(clustertwo[m]-1)]); /* number of measurements */
                        clus[(6*leng)+(clustertwo[m]-1)]=clusterone[m]; /* cluster ID */
                        clus[(7*leng)+(clustertwo[m]-1)]=1;            /* merged? */
                        /* delete that pair & all links to second (=merged) cluster & check blacklist  */
                        gone=clustertwo[m];
                        stay=clusterone[m];
                        mergesize=clusterone.size();
                        n=0;
                        while(n<mergesize){ /* remove merges */
                            if( (clusterone[n]==gone) || (clustertwo[n]==gone) ){
                                clusterone.erase(clusterone.begin()+n);
                                clustertwo.erase(clustertwo.begin()+n);
                                clusterdiff.erase(clusterdiff.begin()+n);
                                mergesize=clusterone.size();
                            }else{
                                n++;
                            }
                        }
                        if((blacksize>0) & (mergesize>0)){
                            clusterase.clear(); /* check blacklist-entries */
                            for(k=0;k<blacksize;k++){
                                if(blackone[k]==gone){
                                    clusterase.push_back(blacktwo[k]);
                                }
                                 if(blacktwo[k]==gone){
                                    clusterase.push_back(blackone[k]);
                                }
                            }
                            if(clusterase.size()>0){ /* remove indirect blacklisted links */
                                for(m=0;(unsigned)m<clusterase.size();m++){
                                    mergesize=clusterone.size();
                                    n=0;
                                    while(n<mergesize){
                                        if( (clustertwo[n]==stay) || (clusterone[n]==clusterase[m])){
                                            clusterone.erase(clusterone.begin()+n);
                                            clustertwo.erase(clustertwo.begin()+n);
                                            clusterdiff.erase(clusterdiff.begin()+n);
                                            mergesize=clusterone.size();
                                        }else{
                                            n++;
                                        }
                                    }
                                }
                            }
                        }
                   }
               }
               /* output ******************************************************/
               delete[] at;
               UNPROTECT(10);
               return(clusters);
            }else{
               /* output ******************************************************/
               delete[] at;
               UNPROTECT(10);
               return(clusters);
           }
       }


      /************************************************************************/
      /* density vector *******************************************************/
      /************************************************************************/
      SEXP densvec(     SEXP RT,
                        SEXP mz,
                        SEXP dmz,
                        SEXP ppm,
                        SEXP drt,
                        SEXP scaled
                        ){

           PROTECT(RT = AS_NUMERIC(RT));
           PROTECT(mz = AS_NUMERIC(mz));
           PROTECT(dmz = AS_NUMERIC(dmz));
           PROTECT(ppm = AS_INTEGER(ppm));
           PROTECT(drt = AS_NUMERIC(drt));
           PROTECT(scaled = AS_NUMERIC(scaled));

           int n;
           int ppm2 = INTEGER_VALUE(ppm);
           int leng2 = LENGTH(RT);
           double dmz2 = NUMERIC_VALUE(dmz);
           double drt2 = NUMERIC_VALUE(drt);
           double *ret, *mass;
           ret = NUMERIC_POINTER(RT);
           mass = NUMERIC_POINTER(mz);
           double *scaled2;
           scaled2 = NUMERIC_POINTER(scaled);
           float *densit;
           densit = new float[leng2];
           for(n=0;n<leng2;n++){
               densit[n]=0;
           }

           densvector(leng2, densit, ret, mass, drt2, dmz2, ppm2, scaled2);

           SEXP outit;
           PROTECT(outit = NEW_NUMERIC(leng2));
           double *at;
           at = NUMERIC_POINTER(outit);
           for(n=0;n<leng2;n++){
               *(at+n) = densit[n];
           }
           SETLENGTH(outit, leng2);
           delete[] densit;
           UNPROTECT(7);
           return outit;

      }


      /************************************************************************/
      /* gap filling, mean for same RT ****************************************/
      /************************************************************************/
       SEXP gapfill(  SEXP RT, /* must be sorted */
                      SEXP intens,
                      SEXP intensorder,
                      SEXP mz,
                      SEXP index,
                      SEXP scans,
                      SEXP drtsmall
                      ){

           PROTECT(RT = AS_NUMERIC(RT));
           PROTECT(intens = AS_NUMERIC(intens));
           PROTECT(intensorder = AS_INTEGER(intensorder));
           PROTECT(mz = AS_NUMERIC(mz));
           PROTECT(index = AS_NUMERIC(index));
           PROTECT(scans = AS_NUMERIC(scans));
           PROTECT(drtsmall = AS_NUMERIC(drtsmall));

           double drt = NUMERIC_VALUE(drtsmall);
           int leng2 = LENGTH(RT);
           int leng3 = LENGTH(scans);

           int n,m,l,p,k;
           double meanint;

           double *ret, *inte, *mass, *scanned, *ind;
           ret = NUMERIC_POINTER(RT);
           inte = NUMERIC_POINTER(intens);
           mass = NUMERIC_POINTER(mz);
           scanned = NUMERIC_POINTER(scans);
           ind = NUMERIC_POINTER(index);
           int *ordret;
           ordret = INTEGER_POINTER(intensorder);

           SEXP ans;
           PROTECT(ans = allocMatrix(REALSXP, leng3, 10));
           double *rans;
           rans = REAL(ans);
           for(m=0;m<leng3;m++){
               rans[m]=0; /* mass */
               rans[m+(1*leng3)]=0; /* intensity */
               rans[m+(2*leng3)]=*(scanned+m); /* RT */
               rans[m+(3*leng3)]=0; /* index */
               rans[m+(4*leng3)]=0; /* Filter */
               rans[m+(5*leng3)]=0; /* 1st peak-pick */
               rans[m+(6*leng3)]=0; /* with peak criteria */
               rans[m+(7*leng3)]=0; /* baseline #1 */
               rans[m+(8*leng3)]=0; /* baseline #2*/
               rans[m+(9*leng3)]=0; /* 2nd peak-pick */
           }

           int *traced;
           traced = new int[leng3];
           for(m=0;m<leng3;m++){
               traced[m]=0;
           }
           /* for identical RT, fill with intensity-closest *******************/
           l=0;
           p=0;
           for(n=0;n<leng2;n++){
               for(m=l;m<leng3;m++){
                   if(*(ret+(*(ordret+n)-1))==rans[m+(2*leng3)]){
                       if(rans[m+(1*leng3)]!=0){
                           traced[m]=traced[m]+1;
                           p=1;
                       }
                       rans[m]=*(mass+(*(ordret+n)-1));
                       rans[m+(1*leng3)]=*(inte+(*(ordret+n)-1));
                       rans[m+(3*leng3)]=*(ind+(*(ordret+n)-1));
                       l=m;
                       break;
                   };
               };
           };
           if(p!=0){
               for(m=0;m<leng3;m++){
                   if(traced[m]!=0){
                       l=0;
                       meanint=0;
                       for(n=(m-2);n<=(m+2);n++){
                           if((n>=0)&(n<leng3)){
                               if((traced[n]==0)&(rans[n+(1*leng3)]!=0)){
                                   meanint=meanint+rans[n+(1*leng3)];
                                   l++;
                               }
                           }
                       }
                       if(l!=0){
                           meanint=meanint/l;
                           for(k=0;k<leng2;k++){
                               if(*(ret+k)==rans[m+(2*leng3)]){
                                   if((fabs(*(inte+k)-meanint))<(fabs(rans[m+(1*leng3)]-meanint))){
                                       rans[m]=*(mass+k);
                                       rans[m+(1*leng3)]=*(inte+k);
                                       rans[m+(3*leng3)]=*(ind+k);
                                   }
                               }
                           }
                       }
                   }
               }
           }
           /* gap-filling - interpolation *************************************/
           k=0;
           for(m=k;m<(leng3-2);m++){
               if(rans[m]!=0){
                   for(n=(m+1);n<leng3;n++){
                       if(rans[n]!=0){
                           break;
                       }
                   }
                   if((n-m)>1){ /* at least one measurement in between */
                       if(fabs(rans[n+(2*leng3)]-rans[m+(2*leng3)])<=drt){
                           for(l=(m+1);l<n;l++){
                               rans[l+(1*leng3)]=rans[m+(1*leng3)]+(
                                   ((rans[n+(1*leng3)]-rans[m+(1*leng3)]))/
                                   (fabs(rans[n+(2*leng3)]-rans[m+(2*leng3)]))*
                                   (fabs(rans[l+(2*leng3)]-rans[m+(2*leng3)]))
                               );
                           }
                       }
                       k=n;
                   }else{
                       k=n;
                   }
               }
           }
           delete[] traced;
           UNPROTECT(8);
           return ans;
      }


      /************************************************************************/
      /* peakpicking **********************************************************/
      /************************************************************************/
       SEXP pickpeak(   SEXP out1, /* generated in gapfill */
                        SEXP drtlarge,
                        SEXP drttotal,
                        SEXP minpeak,
                        SEXP recurs,
                        SEXP weight,
                        SEXP SB,
                        SEXP SN,
                        SEXP minint,
                        SEXP upint,
                        SEXP ended,
                        SEXP win
                      ){

           PROTECT(out1 = AS_NUMERIC(out1));
           PROTECT(drtlarge = AS_NUMERIC(drtlarge));
           PROTECT(drttotal = AS_NUMERIC(drttotal));
           PROTECT(minpeak = AS_INTEGER(minpeak));
           PROTECT(recurs = AS_INTEGER(recurs));
           PROTECT(weight = AS_INTEGER(weight));
           PROTECT(SB = AS_NUMERIC(SB));
           PROTECT(SN = AS_NUMERIC(SN));
           PROTECT(minint = AS_NUMERIC(minint));
           PROTECT(upint = AS_NUMERIC(upint));
           PROTECT(ended = AS_INTEGER(ended));
           PROTECT(win = AS_INTEGER(win));

           double drt2 = NUMERIC_VALUE(drtlarge);
           double drt4 = NUMERIC_VALUE(drttotal);
           double weight2 = NUMERIC_VALUE(weight);
           double minint2 = NUMERIC_VALUE(minint);
           double upint2 = NUMERIC_VALUE(upint);
           double SB2 = NUMERIC_VALUE(SB);
           double SN2 = NUMERIC_VALUE(SN);
           int win2 = INTEGER_VALUE(win);
           int ended2 = INTEGER_VALUE(ended);
           int minpeak2 = INTEGER_VALUE(minpeak);
           int recurs2 = INTEGER_VALUE(recurs);
           int leng3 = LENGTH(out1);
           leng3=leng3/10;

           double *rans;
           rans = REAL(out1);

           double often=0;

           /* 1st recursive peak detection ************************************/
           peakdetect(leng3,&often,recurs2,rans,weight2,drt4);
           if(often>0){
               /* check if peak-criteria fulfilled & filter *******************/
               peakcrit1(rans,leng3,minpeak2,SB2,minint2,upint2,ended2,drt2,&often);
           };
           if(often>0){
               /* subtract + interpolate + get baseline ***********************/
               /* 2nd peak criteria check *************************************/
               peakcrit2(rans,leng3,minpeak2,minint2,upint2,win2,often,SN2);
           }
           UNPROTECT(12);
           return out1;

      }


      /************************************************************************/
      /* generate final peaklist **********************************************/
      /************************************************************************/
      SEXP picklist(  SEXP out2,
                      SEXP peaknumb,
                      SEXP EICnumb
                      ){

           PROTECT(out2 = AS_NUMERIC(out2));
           PROTECT(peaknumb = AS_INTEGER(peaknumb));
           PROTECT(EICnumb = AS_INTEGER(EICnumb));

           int peaknumb2 = INTEGER_VALUE(peaknumb);
           int EICnumb2 = INTEGER_VALUE(EICnumb);
           double *rans;
           rans = REAL(out2);
           int leng3 = LENGTH(out2);
           leng3=leng3/10;
           int n,m,p,from,to;
           double meanmass,varmass,maxint,sumint,minRT,maxRT,intRT,numpeak;
           double leng4=0;
           for(m=0;m<leng3;m++){
               if(leng4<rans[m+(9*leng3)]){
                   leng4=rans[m+(9*leng3)];
               }
           }
           SEXP out3;
           PROTECT(out3 = allocMatrix(REALSXP, int(leng4), 10));
           double *rans2;
           rans2 = REAL(out3);

           for(m=1;m<=int(leng4);m++){
                      /* which peaks in that group, from->to */
                      from=0;
                      to=0;
                      for(n=0;n<leng3;n++){
                          if(rans[n+(9*leng3)]==double(m)){
                              if(from==0){from=n;};
                              to=n;
                          }else{
                              if(from!=0){break;};
                          };
                      };
                      /* extract peak data ************************************/
                      meanmass=0;varmass=0;maxint=0;sumint=0;minRT=0;maxRT=0;intRT=0;
                      numpeak=0;
                      for(p=from;p<=to;p++){
                          if(rans[p]!=0){
                              meanmass=(meanmass+rans[p]);
                              numpeak++;
                          }
                          sumint=(sumint+rans[p+(8*leng3)]);
                          if(rans[p+(8*leng3)]>maxint){
                              maxint=rans[p+(8*leng3)];
                              intRT=rans[p+(2*leng3)];
                          };
                      };
                      meanmass=(meanmass/numpeak);
                      minRT=rans[from+(2*leng3)];
                      maxRT=rans[to+(2*leng3)];
                      for(p=from;p<=to;p++){
                          if(rans[p]!=0){
                              varmass=(varmass+pow((meanmass-rans[p]),2));
                          }
                      };
                      varmass=(varmass/numpeak);
                      rans2[(m-1)+(0*int(leng4))]=meanmass;
                      rans2[(m-1)+(1*int(leng4))]=varmass;
                      rans2[(m-1)+(2*int(leng4))]=maxint;
                      rans2[(m-1)+(3*int(leng4))]=sumint;
                      rans2[(m-1)+(4*int(leng4))]=intRT;
                      rans2[(m-1)+(5*int(leng4))]=minRT;
                      rans2[(m-1)+(6*int(leng4))]=maxRT;
                      rans2[(m-1)+(7*int(leng4))]=double(m+peaknumb2);
                      rans2[(m-1)+(8*int(leng4))]=double(EICnumb2);
                      rans2[(m-1)+(9*int(leng4))]=0;
           };

           UNPROTECT(4);
           return(out3);
      }


      /************************************************************************/
      /* extract measurements for plotting ************************************/
      /************************************************************************/
      SEXP plotit(      SEXP RTlim_low,
                        SEXP RTlim_up,
                        SEXP mzlim_low,
                        SEXP mzlim_up,
                        SEXP mz,
                        SEXP RT,
                        SEXP intensity,
                        SEXP color1,
                        SEXP color2,
                        SEXP color3,
                        SEXP whatcolor
                        ){

           PROTECT(RTlim_low = AS_NUMERIC(RTlim_low));
           PROTECT(RTlim_up = AS_NUMERIC(RTlim_up));
           PROTECT(mzlim_low = AS_NUMERIC(mzlim_low));
           PROTECT(mzlim_up = AS_NUMERIC(mzlim_up));
           PROTECT(RT = AS_NUMERIC(RT));
           PROTECT(mz = AS_NUMERIC(mz));
           PROTECT(intensity = AS_NUMERIC(intensity));
           PROTECT(color1 = AS_NUMERIC(color1));
           PROTECT(color2 = AS_NUMERIC(color2));
           PROTECT(color3 = AS_NUMERIC(color3));
           PROTECT(whatcolor = AS_INTEGER(whatcolor));
           double RTlim_low2 = NUMERIC_VALUE(RTlim_low);
           double RTlim_up2 = NUMERIC_VALUE(RTlim_up);
           double mzlim_low2 = NUMERIC_VALUE(mzlim_low);
           double mzlim_up2 = NUMERIC_VALUE(mzlim_up);
           double whatcol = INTEGER_VALUE(whatcolor);
           double *ret, *mass, *inte, maxint=0;
           ret = NUMERIC_POINTER(RT);
           mass = NUMERIC_POINTER(mz);
           inte = NUMERIC_POINTER(intensity);
           double *colorit1,*colorit2,*colorit3;
           colorit1 = NUMERIC_POINTER(color1);
           colorit2 = NUMERIC_POINTER(color2);
           colorit3 = NUMERIC_POINTER(color3);
           int leng2 = LENGTH(RT);
           int n,counter1=0,counter2=0;
           int *does;
           does = new int[leng2];

           for(n=0;n<leng2;n++){
               if((*(mass+n)>=mzlim_low2) & (*(mass+n)<=mzlim_up2)){
                   if((*(ret+n)>=RTlim_low2) & (*(ret+n)<=RTlim_up2)){
                       does[n]=1;
                       if(*(inte+n)>maxint){
                           maxint=*(inte+n);
                       }
                       counter1++;
                   }else{does[n]=0;};
               }else{does[n]=0;};
           }

           if(counter1>0){
               SEXP ans;
               PROTECT(ans = allocMatrix(REALSXP, counter1, 5));
               double *rans;
               rans = REAL(ans);
               if(whatcol==1){
                   for(n=0;n<leng2;n++){
                       if(does[n]==1){
                           rans[counter2]=*(mass+n);
                           rans[counter2+counter1]=*(inte+n);
                           rans[counter2+(2*counter1)]=*(ret+n);
                           rans[counter2+(3*counter1)]=*(colorit1+n)+1;
                           rans[counter2+(4*counter1)]=0;
                           counter2++;
                       }
                   }
               }
               if(whatcol==2){
                   for(n=0;n<leng2;n++){
                       if(does[n]==1){
                           rans[counter2]=*(mass+n);
                           rans[counter2+counter1]=*(inte+n);
                           rans[counter2+(2*counter1)]=*(ret+n);
                           rans[counter2+(3*counter1)]=*(colorit2+n)+1;
                           rans[counter2+(4*counter1)]=0;
                           counter2++;
                       }
                   }
               }
               if(whatcol==3){
                   for(n=0;n<leng2;n++){
                       if(does[n]==1){
                           rans[counter2]=*(mass+n);
                           rans[counter2+counter1]=*(inte+n);
                           rans[counter2+(2*counter1)]=*(ret+n);
                           rans[counter2+(3*counter1)]=*(colorit3+n)+1;
                           rans[counter2+(4*counter1)]=0;
                           counter2++;
                       }
                   }
               }
               if(whatcol==4){
                   for(n=0;n<leng2;n++){
                       if(does[n]==1){
                           rans[counter2]=*(mass+n);
                           rans[counter2+counter1]=*(inte+n);
                           rans[counter2+(2*counter1)]=*(ret+n);
                           rans[counter2+(3*counter1)]=1;
                           rans[counter2+(4*counter1)]=1;
                           counter2++;
                       }
                   }
               }
               rans[(4*counter1)]=maxint;
               rans[(4*counter1)+1]=counter1;
			   delete does;
               SETLENGTH(ans, counter1*5);
               UNPROTECT(12);
               return ans;
           }else{
           	   delete does;
               UNPROTECT(11);
               return R_NilValue;
           }

      }


      /************************************************************************/
      /* bin measurements for plotting - RT ***********************************/
      /************************************************************************/
      SEXP binRT(  SEXP RT,
                   SEXP intensity,
                   SEXP scantimes,
                   SEXP colorit,
                   SEXP what
                   ){

           PROTECT(RT = AS_NUMERIC(RT));
           PROTECT(intensity = AS_NUMERIC(intensity));
           PROTECT(scantimes = AS_NUMERIC(scantimes));
           PROTECT(colorit = AS_NUMERIC(colorit));
           PROTECT(what = AS_INTEGER(what));
           int what2 = INTEGER_VALUE(what);
           int leng2 = LENGTH(RT);
           int leng3 = LENGTH(scantimes);
           double *ret, *intens, *scans, *col;
           ret = NUMERIC_POINTER(RT);
           intens = NUMERIC_POINTER(intensity);
           scans = NUMERIC_POINTER(scantimes);
           col = NUMERIC_POINTER(colorit);
           int n,m,l,counter;
           SEXP ans;
           PROTECT(ans = allocMatrix(REALSXP, leng3, 1));
           double *rans;
           rans = REAL(ans);
           for(m=0;m<leng3;m++){
               rans[m]=0;
           };

           /* find entries, set to most intensive */
           if(what2==1){
               l=0;
               counter=0;
               for(n=0;n<leng2;n++){ /* over data*/
                   if(*(col+n)>1){
                       for(m=l;m<leng3;m++){ /* over scans */
                           if(*(ret+n)==*(scans+m)){
                               if(rans[m]==0){
                                    counter++;
                               }
                               rans[m]=rans[m]+*(intens+n);
                               l=m;
                               break;
                           }
                       }
                   }
               }
           }else{
               l=0;
               counter=0;
               for(n=0;n<leng2;n++){ /* over data*/
                   if(*(col+n)>1){
                       for(m=l;m<leng3;m++){ /* over scans */
                           if(*(ret+n)==*(scans+m)){
                               if(rans[m]==0){
                                    counter++;
                               }
                               if(rans[m]<*(intens+n)){
                                   rans[m]=*(intens+n);
                               }
                               l=m;
                               break;
                           }
                       }
                   }
               }
           }

           /* omit 0-entries in ans2, i.e., resize to within RTlimit set in R */
           SEXP ans2;
           PROTECT(ans2 = allocMatrix(REALSXP, counter, 2));
           double *rans2;
           rans2 = REAL(ans2);
           n=0;
           for(m=0;m<leng3;m++){
               if(rans[m]!=0){
                   rans2[n]=rans[m];
                   rans2[counter+n]=*(scans+m);
                   n++;
               }
           }
           UNPROTECT(7);
           return ans2;
      }


      /************************************************************************/
      /* bin measurements for plotting - mass *********************************/
      /************************************************************************/
      SEXP binmz(  SEXP mz,
                   SEXP intensity,
                   SEXP binmzs,
                   SEXP colorit
                   ){

           PROTECT(mz = AS_NUMERIC(mz));
           PROTECT(intensity = AS_NUMERIC(intensity));
           PROTECT(binmzs = AS_NUMERIC(binmzs));
           PROTECT(colorit = AS_NUMERIC(colorit));
           int leng2 = LENGTH(mz);
           int leng3 = LENGTH(binmzs);
           double *mass, *intens, *binit, *col;
           mass = NUMERIC_POINTER(mz);
           intens = NUMERIC_POINTER(intensity);
           binit = NUMERIC_POINTER(binmzs);
           col = NUMERIC_POINTER(colorit);
           int n,m,l,counter;
           SEXP ans;
           PROTECT(ans = allocMatrix(REALSXP, leng3-1, 2));
           double *rans;
           rans = REAL(ans);
           for(m=0;m<(leng3-1);m++){
               rans[m]=0;
               rans[m+leng3-1]=0;
           };

           l=0;
           counter=0;
           for(n=0;n<leng2;n++){
               if(*(col+n)>1){
                   for(m=l;m<(leng3-1);m++){
                       if((*(mass+n)>=*(binit+m))&(*(mass+n)<*(binit+m+1))){
                           if(rans[m]==0){counter++;};
                               if(rans[m]<*(intens+n)){
                                   rans[m]=*(intens+n);
                                   rans[m+leng3-1]=((*(binit+m)+*(binit+m+1))/2);
                               }
                           l=m;
                           break;
                       }
                   }
               }
           }
           /* omit 0-entries in ans2, i.e., resize to within RTlimit set in R */
           SEXP ans2;
           PROTECT(ans2 = allocMatrix(REALSXP, counter, 2));
           double *rans2;
           rans2 = REAL(ans2);
           n=0;
           for(m=0;m<(leng3-1);m++){
               if(rans[m]!=0){
                   rans2[n]=rans[m];
                   rans2[counter+n]=rans[m+leng3-1];
                   n++;
               }
           }
           UNPROTECT(6);
           return ans2;
     }


}




