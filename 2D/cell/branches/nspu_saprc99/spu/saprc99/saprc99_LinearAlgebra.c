/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Linear Algebra Data and Routines File                            */
/*                                                                  */
/* Generated by KPP-2.2 symbolic chemistry Kinetics PreProcessor    */
/*       (http://www.cs.vt.edu/~asandu/Software/KPP)                */
/* KPP is distributed under GPL, the general public licence         */
/*       (http://www.gnu.org/copyleft/gpl.html)                     */
/* (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa           */
/* (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech            */
/*     With important contributions from:                           */
/*        M. Damian, Villanova University, USA                      */
/*        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany */
/*                                                                  */
/* File                 : saprc99_LinearAlgebra.c                   */
/* Time                 : Wed Jun 11 10:09:48 2008                  */
/* Working directory    : /Users/jlinford/workspace/saprc99         */
/* Equation file        : saprc99.kpp                               */
/* Output root filename : saprc99                                   */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#define USE_SDK_BLAS    0

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "saprc99_Parameters.h"
#include "saprc99_Global.h"
#include "saprc99_Sparse.h"

/* For SPU-optimized BLAS */
#if USE_SDK_BLAS == 1
#include <blas.h>
#endif


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* SPARSE_UTIL - SPARSE utility functions                           */
/*   Arguments :                                                    */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


int KppDecomp( double *JVS )
{
double W[74];
double a;
int k, kk, j, jj;
    
  for( k = 0; k < 74; k++ ) {
    if( JVS[ LU_DIAG[k] ] == 0.0 ) return k+1;
    for( kk = LU_CROW[k]; kk < LU_CROW[k+1]; kk++ )
      W[ LU_ICOL[kk] ] = JVS[kk];
    for( kk = LU_CROW[k]; kk < LU_DIAG[k]; kk++ ) {
      j = LU_ICOL[kk];
      a = -W[j] / JVS[ LU_DIAG[j] ];
      W[j] = -a;
      for( jj = LU_DIAG[j]+1; jj < LU_CROW[j+1]; jj++ )
        W[ LU_ICOL[jj] ] += a*JVS[jj];
    }
    for( kk = LU_CROW[k]; kk < LU_CROW[k+1]; kk++ )
      JVS[kk] = W[ LU_ICOL[kk] ];
  }
  return 0;
}

/* End of SPARSE_UTIL function                                      */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* KppSolve - sparse back substitution                              */
/*   Arguments :                                                    */
/*      JVS       - sparse Jacobian of variables                    */
/*      X         - Vector for variables                            */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

void KppSolve( 
  double JVS[],                           /* sparse Jacobian of variables */
  double X[]                              /* Vector for variables */
)
{
  X[21] = X[21]-JVS[91]*X[20];
  X[31] = X[31]-JVS[125]*X[23]-JVS[126]*X[30];
  X[32] = X[32]-JVS[129]*X[23]-JVS[130]*X[30];
  X[33] = X[33]-JVS[133]*X[23]-JVS[134]*X[30];
  X[34] = X[34]-JVS[138]*X[23]-JVS[139]*X[30];
  X[35] = X[35]-JVS[143]*X[27];
  X[37] = X[37]-JVS[154]*X[23]-JVS[155]*X[30];
  X[38] = X[38]-JVS[161]*X[30];
  X[39] = X[39]-JVS[167]*X[19]-JVS[168]*X[29]-JVS[169]*X[31]-JVS[170]
         *X[32]-JVS[171]*X[34];
  X[40] = X[40]-JVS[190]*X[23]-JVS[191]*X[30]-JVS[192]*X[31]-JVS[193]
         *X[32]-JVS[194]*X[33];
  X[41] = X[41]-JVS[202]*X[19]-JVS[203]*X[20]-JVS[204]*X[21]-JVS[205]
         *X[22]-JVS[206]*X[29];
  X[42] = X[42]-JVS[216]*X[17]-JVS[217]*X[33]-JVS[218]*X[35]-JVS[219]
         *X[37]-JVS[220]*X[38]-JVS[221]*X[40];
  X[44] = X[44]-JVS[242]*X[19]-JVS[243]*X[23]-JVS[244]*X[30]-JVS[245]
         *X[31]-JVS[246]*X[32]-JVS[247]*X[34]-JVS[248]*X[38]-JVS[249]
         *X[43];
  X[45] = X[45]-JVS[259]*X[33]-JVS[260]*X[38];
  X[47] = X[47]-JVS[276]*X[20]-JVS[277]*X[22]-JVS[278]*X[29]-JVS[279]
         *X[31]-JVS[280]*X[32]-JVS[281]*X[41]-JVS[282]*X[46];
  X[49] = X[49]-JVS[310]*X[46];
  X[51] = X[51]-JVS[322]*X[46];
  X[53] = X[53]-JVS[334]*X[46]-JVS[335]*X[52];
  X[54] = X[54]-JVS[341]*X[10]-JVS[342]*X[20]-JVS[343]*X[22]-JVS[344]
         *X[29]-JVS[345]*X[43]-JVS[346]*X[50]-JVS[347]*X[51]-JVS[348]
         *X[52];
  X[55] = X[55]-JVS[363]*X[19]-JVS[364]*X[20]-JVS[365]*X[22]-JVS[366]
         *X[25]-JVS[367]*X[26]-JVS[368]*X[28]-JVS[369]*X[29]-JVS[370]
         *X[41]-JVS[371]*X[43]-JVS[372]*X[44]-JVS[373]*X[45]-JVS[374]
         *X[46]-JVS[375]*X[48]-JVS[376]*X[49]-JVS[377]*X[50]-JVS[378]
         *X[51]-JVS[379]*X[52]-JVS[380]*X[53];
  X[56] = X[56]-JVS[399]*X[21]-JVS[400]*X[48]-JVS[401]*X[50]-JVS[402]
         *X[51]-JVS[403]*X[52];
  X[57] = X[57]-JVS[412]*X[9]-JVS[413]*X[43]-JVS[414]*X[46]-JVS[415]
         *X[48]-JVS[416]*X[49]-JVS[417]*X[50]-JVS[418]*X[52]-JVS[419]
         *X[53];
  X[58] = X[58]-JVS[426]*X[19]-JVS[427]*X[20]-JVS[428]*X[22]-JVS[429]
         *X[29]-JVS[430]*X[31]-JVS[431]*X[32]-JVS[432]*X[34]-JVS[433]
         *X[36]-JVS[434]*X[43]-JVS[435]*X[48]-JVS[436]*X[49]-JVS[437]
         *X[50]-JVS[438]*X[51]-JVS[439]*X[52]-JVS[440]*X[53]-JVS[441]
         *X[56]-JVS[442]*X[57];
  X[59] = X[59]-JVS[454]*X[20]-JVS[455]*X[22]-JVS[456]*X[29]-JVS[457]
         *X[49]-JVS[458]*X[50]-JVS[459]*X[51]-JVS[460]*X[52]-JVS[461]
         *X[53]-JVS[462]*X[56]-JVS[463]*X[57];
  X[60] = X[60]-JVS[474]*X[22]-JVS[475]*X[29]-JVS[476]*X[30]-JVS[477]
         *X[46]-JVS[478]*X[48]-JVS[479]*X[50]-JVS[480]*X[51]-JVS[481]
         *X[52]-JVS[482]*X[53]-JVS[483]*X[56]-JVS[484]*X[57];
  X[61] = X[61]-JVS[497]*X[34]-JVS[498]*X[43]-JVS[499]*X[46]-JVS[500]
         *X[48]-JVS[501]*X[49]-JVS[502]*X[50]-JVS[503]*X[51]-JVS[504]
         *X[52]-JVS[505]*X[53]-JVS[506]*X[57];
  X[62] = X[62]-JVS[517]*X[8]-JVS[518]*X[16]-JVS[519]*X[18]-JVS[520]
         *X[19]-JVS[521]*X[23]-JVS[522]*X[24]-JVS[523]*X[25]-JVS[524]
         *X[26]-JVS[525]*X[27]-JVS[526]*X[28]-JVS[527]*X[30]-JVS[528]
         *X[31]-JVS[529]*X[32]-JVS[530]*X[34]-JVS[531]*X[35]-JVS[532]
         *X[36]-JVS[533]*X[39]-JVS[534]*X[40]-JVS[535]*X[43]-JVS[536]
         *X[44]-JVS[537]*X[45]-JVS[538]*X[46]-JVS[539]*X[48]-JVS[540]
         *X[49]-JVS[541]*X[50]-JVS[542]*X[51]-JVS[543]*X[52]-JVS[544]
         *X[53]-JVS[545]*X[54]-JVS[546]*X[55]-JVS[547]*X[56]-JVS[548]
         *X[57]-JVS[549]*X[58]-JVS[550]*X[59]-JVS[551]*X[60]-JVS[552]
         *X[61];
  X[63] = X[63]-JVS[565]*X[19]-JVS[566]*X[20]-JVS[567]*X[22]-JVS[568]
         *X[23]-JVS[569]*X[29]-JVS[570]*X[30]-JVS[571]*X[46]-JVS[572]
         *X[48]-JVS[573]*X[50]-JVS[574]*X[51]-JVS[575]*X[52]-JVS[576]
         *X[53]-JVS[577]*X[56]-JVS[578]*X[57]-JVS[579]*X[58]-JVS[580]
         *X[59]-JVS[581]*X[60]-JVS[582]*X[61]-JVS[583]*X[62];
  X[64] = X[64]-JVS[595]*X[15]-JVS[596]*X[46]-JVS[597]*X[49]-JVS[598]
         *X[51]-JVS[599]*X[52]-JVS[600]*X[53]-JVS[601]*X[57]-JVS[602]
         *X[61]-JVS[603]*X[62]-JVS[604]*X[63];
  X[65] = X[65]-JVS[615]*X[21]-JVS[616]*X[25]-JVS[617]*X[29]-JVS[618]
         *X[41]-JVS[619]*X[43]-JVS[620]*X[46]-JVS[621]*X[48]-JVS[622]
         *X[50]-JVS[623]*X[52]-JVS[624]*X[53]-JVS[625]*X[54]-JVS[626]
         *X[56]-JVS[627]*X[57]-JVS[628]*X[58]-JVS[629]*X[59]-JVS[630]
         *X[60]-JVS[631]*X[61]-JVS[632]*X[62]-JVS[633]*X[63]-JVS[634]
         *X[64];
  X[66] = X[66]-JVS[644]*X[14]-JVS[645]*X[37]-JVS[646]*X[52]-JVS[647]
         *X[57]-JVS[648]*X[61]-JVS[649]*X[62]-JVS[650]*X[63]-JVS[651]
         *X[64]-JVS[652]*X[65];
  X[67] = X[67]-JVS[661]*X[10]-JVS[662]*X[19]-JVS[663]*X[20]-JVS[664]
         *X[22]-JVS[665]*X[23]-JVS[666]*X[29]-JVS[667]*X[30]-JVS[668]
         *X[31]-JVS[669]*X[32]-JVS[670]*X[33]-JVS[671]*X[34]-JVS[672]
         *X[36]-JVS[673]*X[38]-JVS[674]*X[43]-JVS[675]*X[45]-JVS[676]
         *X[46]-JVS[677]*X[48]-JVS[678]*X[49]-JVS[679]*X[50]-JVS[680]
         *X[51]-JVS[681]*X[52]-JVS[682]*X[53]-JVS[683]*X[56]-JVS[684]
         *X[57]-JVS[685]*X[58]-JVS[686]*X[59]-JVS[687]*X[60]-JVS[688]
         *X[61]-JVS[689]*X[62]-JVS[690]*X[63]-JVS[691]*X[64]-JVS[692]
         *X[65]-JVS[693]*X[66];
  X[68] = X[68]-JVS[701]*X[18]-JVS[702]*X[26]-JVS[703]*X[47]-JVS[704]
         *X[48]-JVS[705]*X[50]-JVS[706]*X[52]-JVS[707]*X[53]-JVS[708]
         *X[55]-JVS[709]*X[56]-JVS[710]*X[57]-JVS[711]*X[59]-JVS[712]
         *X[60]-JVS[713]*X[61]-JVS[714]*X[62]-JVS[715]*X[63]-JVS[716]
         *X[64]-JVS[717]*X[65]-JVS[718]*X[66]-JVS[719]*X[67];
  X[69] = X[69]-JVS[726]*X[12]-JVS[727]*X[13]-JVS[728]*X[14]-JVS[729]
         *X[15]-JVS[730]*X[17]-JVS[731]*X[18]-JVS[732]*X[21]-JVS[733]
         *X[24]-JVS[734]*X[26]-JVS[735]*X[27]-JVS[736]*X[35]-JVS[737]
         *X[42]-JVS[738]*X[44]-JVS[739]*X[45]-JVS[740]*X[46]-JVS[741]
         *X[47]-JVS[742]*X[48]-JVS[743]*X[49]-JVS[744]*X[50]-JVS[745]
         *X[51]-JVS[746]*X[52]-JVS[747]*X[53]-JVS[748]*X[54]-JVS[749]
         *X[55]-JVS[750]*X[56]-JVS[751]*X[57]-JVS[752]*X[58]-JVS[753]
         *X[59]-JVS[754]*X[60]-JVS[755]*X[61]-JVS[756]*X[62]-JVS[757]
         *X[63]-JVS[758]*X[64]-JVS[759]*X[65]-JVS[760]*X[66]-JVS[761]
         *X[67]-JVS[762]*X[68];
  X[70] = X[70]-JVS[768]*X[17]-JVS[769]*X[24]-JVS[770]*X[33]-JVS[771]
         *X[35]-JVS[772]*X[37]-JVS[773]*X[38]-JVS[774]*X[40]-JVS[775]
         *X[42]-JVS[776]*X[43]-JVS[777]*X[44]-JVS[778]*X[45]-JVS[779]
         *X[46]-JVS[780]*X[47]-JVS[781]*X[48]-JVS[782]*X[49]-JVS[783]
         *X[50]-JVS[784]*X[51]-JVS[785]*X[52]-JVS[786]*X[53]-JVS[787]
         *X[54]-JVS[788]*X[55]-JVS[789]*X[56]-JVS[790]*X[57]-JVS[791]
         *X[58]-JVS[792]*X[59]-JVS[793]*X[60]-JVS[794]*X[61]-JVS[795]
         *X[62]-JVS[796]*X[63]-JVS[797]*X[64]-JVS[798]*X[65]-JVS[799]
         *X[66]-JVS[800]*X[67]-JVS[801]*X[68]-JVS[802]*X[69];
  X[71] = X[71]-JVS[807]*X[11]-JVS[808]*X[12]-JVS[809]*X[23]-JVS[810]
         *X[29]-JVS[811]*X[31]-JVS[812]*X[32]-JVS[813]*X[40]-JVS[814]
         *X[41]-JVS[815]*X[48]-JVS[816]*X[49]-JVS[817]*X[50]-JVS[818]
         *X[51]-JVS[819]*X[52]-JVS[820]*X[53]-JVS[821]*X[54]-JVS[822]
         *X[56]-JVS[823]*X[57]-JVS[824]*X[58]-JVS[825]*X[59]-JVS[826]
         *X[60]-JVS[827]*X[61]-JVS[828]*X[62]-JVS[829]*X[63]-JVS[830]
         *X[64]-JVS[831]*X[65]-JVS[832]*X[66]-JVS[833]*X[67]-JVS[834]
         *X[68]-JVS[835]*X[69]-JVS[836]*X[70];
  X[72] = X[72]-JVS[840]*X[13]-JVS[841]*X[44]-JVS[842]*X[45]-JVS[843]
         *X[48]-JVS[844]*X[49]-JVS[845]*X[51]-JVS[846]*X[52]-JVS[847]
         *X[53]-JVS[848]*X[57]-JVS[849]*X[58]-JVS[850]*X[59]-JVS[851]
         *X[60]-JVS[852]*X[61]-JVS[853]*X[62]-JVS[854]*X[63]-JVS[855]
         *X[64]-JVS[856]*X[65]-JVS[857]*X[66]-JVS[858]*X[67]-JVS[859]
         *X[68]-JVS[860]*X[69]-JVS[861]*X[70]-JVS[862]*X[71];
  X[73] = X[73]-JVS[865]*X[8]-JVS[866]*X[9]-JVS[867]*X[10]-JVS[868]
         *X[16]-JVS[869]*X[18]-JVS[870]*X[19]-JVS[871]*X[20]-JVS[872]
         *X[22]-JVS[873]*X[23]-JVS[874]*X[24]-JVS[875]*X[25]-JVS[876]
         *X[28]-JVS[877]*X[29]-JVS[878]*X[30]-JVS[879]*X[31]-JVS[880]
         *X[32]-JVS[881]*X[33]-JVS[882]*X[34]-JVS[883]*X[36]-JVS[884]
         *X[37]-JVS[885]*X[38]-JVS[886]*X[39]-JVS[887]*X[40]-JVS[888]
         *X[41]-JVS[889]*X[42]-JVS[890]*X[43]-JVS[891]*X[44]-JVS[892]
         *X[45]-JVS[893]*X[46]-JVS[894]*X[48]-JVS[895]*X[49]-JVS[896]
         *X[50]-JVS[897]*X[51]-JVS[898]*X[52]-JVS[899]*X[53]-JVS[900]
         *X[54]-JVS[901]*X[55]-JVS[902]*X[56]-JVS[903]*X[57]-JVS[904]
         *X[58]-JVS[905]*X[59]-JVS[906]*X[60]-JVS[907]*X[61]-JVS[908]
         *X[62]-JVS[909]*X[63]-JVS[910]*X[64]-JVS[911]*X[65]-JVS[912]
         *X[66]-JVS[913]*X[67]-JVS[914]*X[68]-JVS[915]*X[69]-JVS[916]
         *X[70]-JVS[917]*X[71]-JVS[918]*X[72];
  X[73] = X[73]/JVS[919];
  X[72] = (X[72]-JVS[864]*X[73])/(JVS[863]);
  X[71] = (X[71]-JVS[838]*X[72]-JVS[839]*X[73])/(JVS[837]);
  X[70] = (X[70]-JVS[804]*X[71]-JVS[805]*X[72]-JVS[806]*X[73])
         /(JVS[803]);
  X[69] = (X[69]-JVS[764]*X[70]-JVS[765]*X[71]-JVS[766]*X[72]-JVS[767]
         *X[73])/(JVS[763]);
  X[68] = (X[68]-JVS[721]*X[69]-JVS[722]*X[70]-JVS[723]*X[71]-JVS[724]
         *X[72]-JVS[725]*X[73])/(JVS[720]);
  X[67] = (X[67]-JVS[695]*X[68]-JVS[696]*X[69]-JVS[697]*X[70]-JVS[698]
         *X[71]-JVS[699]*X[72]-JVS[700]*X[73])/(JVS[694]);
  X[66] = (X[66]-JVS[654]*X[67]-JVS[655]*X[68]-JVS[656]*X[69]-JVS[657]
         *X[70]-JVS[658]*X[71]-JVS[659]*X[72]-JVS[660]*X[73])
         /(JVS[653]);
  X[65] = (X[65]-JVS[636]*X[66]-JVS[637]*X[67]-JVS[638]*X[68]-JVS[639]
         *X[69]-JVS[640]*X[70]-JVS[641]*X[71]-JVS[642]*X[72]-JVS[643]
         *X[73])/(JVS[635]);
  X[64] = (X[64]-JVS[606]*X[65]-JVS[607]*X[66]-JVS[608]*X[67]-JVS[609]
         *X[68]-JVS[610]*X[69]-JVS[611]*X[70]-JVS[612]*X[71]-JVS[613]
         *X[72]-JVS[614]*X[73])/(JVS[605]);
  X[63] = (X[63]-JVS[585]*X[64]-JVS[586]*X[65]-JVS[587]*X[66]-JVS[588]
         *X[67]-JVS[589]*X[68]-JVS[590]*X[69]-JVS[591]*X[70]-JVS[592]
         *X[71]-JVS[593]*X[72]-JVS[594]*X[73])/(JVS[584]);
  X[62] = (X[62]-JVS[554]*X[63]-JVS[555]*X[64]-JVS[556]*X[65]-JVS[557]
         *X[66]-JVS[558]*X[67]-JVS[559]*X[68]-JVS[560]*X[69]-JVS[561]
         *X[70]-JVS[562]*X[71]-JVS[563]*X[72]-JVS[564]*X[73])
         /(JVS[553]);
  X[61] = (X[61]-JVS[508]*X[62]-JVS[509]*X[64]-JVS[510]*X[66]-JVS[511]
         *X[68]-JVS[512]*X[69]-JVS[513]*X[70]-JVS[514]*X[71]-JVS[515]
         *X[72]-JVS[516]*X[73])/(JVS[507]);
  X[60] = (X[60]-JVS[486]*X[61]-JVS[487]*X[63]-JVS[488]*X[65]-JVS[489]
         *X[66]-JVS[490]*X[67]-JVS[491]*X[68]-JVS[492]*X[69]-JVS[493]
         *X[70]-JVS[494]*X[71]-JVS[495]*X[72]-JVS[496]*X[73])
         /(JVS[485]);
  X[59] = (X[59]-JVS[465]*X[60]-JVS[466]*X[61]-JVS[467]*X[63]-JVS[468]
         *X[65]-JVS[469]*X[67]-JVS[470]*X[68]-JVS[471]*X[69]-JVS[472]
         *X[70]-JVS[473]*X[73])/(JVS[464]);
  X[58] = (X[58]-JVS[444]*X[59]-JVS[445]*X[60]-JVS[446]*X[61]-JVS[447]
         *X[62]-JVS[448]*X[63]-JVS[449]*X[67]-JVS[450]*X[68]-JVS[451]
         *X[69]-JVS[452]*X[70]-JVS[453]*X[73])/(JVS[443]);
  X[57] = (X[57]-JVS[421]*X[61]-JVS[422]*X[68]-JVS[423]*X[69]-JVS[424]
         *X[70]-JVS[425]*X[73])/(JVS[420]);
  X[56] = (X[56]-JVS[405]*X[57]-JVS[406]*X[61]-JVS[407]*X[63]-JVS[408]
         *X[68]-JVS[409]*X[69]-JVS[410]*X[70]-JVS[411]*X[73])
         /(JVS[404]);
  X[55] = (X[55]-JVS[382]*X[56]-JVS[383]*X[57]-JVS[384]*X[59]-JVS[385]
         *X[60]-JVS[386]*X[61]-JVS[387]*X[62]-JVS[388]*X[63]-JVS[389]
         *X[64]-JVS[390]*X[65]-JVS[391]*X[66]-JVS[392]*X[67]-JVS[393]
         *X[68]-JVS[394]*X[69]-JVS[395]*X[70]-JVS[396]*X[71]-JVS[397]
         *X[72]-JVS[398]*X[73])/(JVS[381]);
  X[54] = (X[54]-JVS[350]*X[56]-JVS[351]*X[57]-JVS[352]*X[58]-JVS[353]
         *X[59]-JVS[354]*X[60]-JVS[355]*X[61]-JVS[356]*X[64]-JVS[357]
         *X[66]-JVS[358]*X[68]-JVS[359]*X[70]-JVS[360]*X[71]-JVS[361]
         *X[72]-JVS[362]*X[73])/(JVS[349]);
  X[53] = (X[53]-JVS[337]*X[57]-JVS[338]*X[61]-JVS[339]*X[70]-JVS[340]
         *X[73])/(JVS[336]);
  X[52] = (X[52]-JVS[330]*X[57]-JVS[331]*X[61]-JVS[332]*X[70]-JVS[333]
         *X[73])/(JVS[329]);
  X[51] = (X[51]-JVS[324]*X[52]-JVS[325]*X[57]-JVS[326]*X[61]-JVS[327]
         *X[70]-JVS[328]*X[73])/(JVS[323]);
  X[50] = (X[50]-JVS[318]*X[57]-JVS[319]*X[61]-JVS[320]*X[70]-JVS[321]
         *X[73])/(JVS[317]);
  X[49] = (X[49]-JVS[312]*X[52]-JVS[313]*X[57]-JVS[314]*X[61]-JVS[315]
         *X[70]-JVS[316]*X[73])/(JVS[311]);
  X[48] = (X[48]-JVS[306]*X[57]-JVS[307]*X[61]-JVS[308]*X[70]-JVS[309]
         *X[73])/(JVS[305]);
  X[47] = (X[47]-JVS[284]*X[48]-JVS[285]*X[50]-JVS[286]*X[52]-JVS[287]
         *X[53]-JVS[288]*X[56]-JVS[289]*X[57]-JVS[290]*X[59]-JVS[291]
         *X[60]-JVS[292]*X[61]-JVS[293]*X[62]-JVS[294]*X[63]-JVS[295]
         *X[64]-JVS[296]*X[65]-JVS[297]*X[66]-JVS[298]*X[67]-JVS[299]
         *X[68]-JVS[300]*X[69]-JVS[301]*X[70]-JVS[302]*X[71]-JVS[303]
         *X[72]-JVS[304]*X[73])/(JVS[283]);
  X[46] = (X[46]-JVS[272]*X[57]-JVS[273]*X[61]-JVS[274]*X[70]-JVS[275]
         *X[73])/(JVS[271]);
  X[45] = (X[45]-JVS[262]*X[62]-JVS[263]*X[64]-JVS[264]*X[66]-JVS[265]
         *X[68]-JVS[266]*X[69]-JVS[267]*X[70]-JVS[268]*X[71]-JVS[269]
         *X[72]-JVS[270]*X[73])/(JVS[261]);
  X[44] = (X[44]-JVS[251]*X[45]-JVS[252]*X[48]-JVS[253]*X[51]-JVS[254]
         *X[57]-JVS[255]*X[61]-JVS[256]*X[62]-JVS[257]*X[70]-JVS[258]
         *X[73])/(JVS[250]);
  X[43] = (X[43]-JVS[238]*X[57]-JVS[239]*X[61]-JVS[240]*X[70]-JVS[241]
         *X[73])/(JVS[237]);
  X[42] = (X[42]-JVS[223]*X[44]-JVS[224]*X[45]-JVS[225]*X[49]-JVS[226]
         *X[51]-JVS[227]*X[52]-JVS[228]*X[53]-JVS[229]*X[54]-JVS[230]
         *X[55]-JVS[231]*X[58]-JVS[232]*X[61]-JVS[233]*X[62]-JVS[234]
         *X[69]-JVS[235]*X[70]-JVS[236]*X[73])/(JVS[222]);
  X[41] = (X[41]-JVS[208]*X[48]-JVS[209]*X[50]-JVS[210]*X[52]-JVS[211]
         *X[56]-JVS[212]*X[61]-JVS[213]*X[69]-JVS[214]*X[70]-JVS[215]
         *X[73])/(JVS[207]);
  X[40] = (X[40]-JVS[196]*X[49]-JVS[197]*X[51]-JVS[198]*X[53]-JVS[199]
         *X[61]-JVS[200]*X[70]-JVS[201]*X[73])/(JVS[195]);
  X[39] = (X[39]-JVS[173]*X[40]-JVS[174]*X[43]-JVS[175]*X[44]-JVS[176]
         *X[46]-JVS[177]*X[48]-JVS[178]*X[49]-JVS[179]*X[50]-JVS[180]
         *X[51]-JVS[181]*X[52]-JVS[182]*X[53]-JVS[183]*X[54]-JVS[184]
         *X[55]-JVS[185]*X[57]-JVS[186]*X[58]-JVS[187]*X[61]-JVS[188]
         *X[70]-JVS[189]*X[73])/(JVS[172]);
  X[38] = (X[38]-JVS[163]*X[45]-JVS[164]*X[62]-JVS[165]*X[70]-JVS[166]
         *X[73])/(JVS[162]);
  X[37] = (X[37]-JVS[157]*X[52]-JVS[158]*X[61]-JVS[159]*X[70]-JVS[160]
         *X[73])/(JVS[156]);
  X[36] = (X[36]-JVS[150]*X[62]-JVS[151]*X[63]-JVS[152]*X[67]-JVS[153]
         *X[73])/(JVS[149]);
  X[35] = (X[35]-JVS[145]*X[45]-JVS[146]*X[62]-JVS[147]*X[69]-JVS[148]
         *X[70])/(JVS[144]);
  X[34] = (X[34]-JVS[141]*X[61]-JVS[142]*X[73])/(JVS[140]);
  X[33] = (X[33]-JVS[136]*X[70]-JVS[137]*X[73])/(JVS[135]);
  X[32] = (X[32]-JVS[132]*X[73])/(JVS[131]);
  X[31] = (X[31]-JVS[128]*X[73])/(JVS[127]);
  X[30] = (X[30]-JVS[124]*X[73])/(JVS[123]);
  X[29] = (X[29]-JVS[122]*X[73])/(JVS[121]);
  X[28] = (X[28]-JVS[117]*X[63]-JVS[118]*X[65]-JVS[119]*X[67]-JVS[120]
         *X[73])/(JVS[116]);
  X[27] = (X[27]-JVS[112]*X[35]-JVS[113]*X[62]-JVS[114]*X[69]-JVS[115]
         *X[70])/(JVS[111]);
  X[26] = (X[26]-JVS[108]*X[55]-JVS[109]*X[62]-JVS[110]*X[68])
         /(JVS[107]);
  X[25] = (X[25]-JVS[104]*X[62]-JVS[105]*X[65]-JVS[106]*X[73])
         /(JVS[103]);
  X[24] = (X[24]-JVS[100]*X[62]-JVS[101]*X[69]-JVS[102]*X[73])
         /(JVS[99]);
  X[23] = (X[23]-JVS[98]*X[73])/(JVS[97]);
  X[22] = (X[22]-JVS[96]*X[73])/(JVS[95]);
  X[21] = (X[21]-JVS[93]*X[69]-JVS[94]*X[73])/(JVS[92]);
  X[20] = (X[20]-JVS[90]*X[73])/(JVS[89]);
  X[19] = (X[19]-JVS[88]*X[73])/(JVS[87]);
  X[18] = (X[18]-JVS[85]*X[68]-JVS[86]*X[73])/(JVS[84]);
  X[17] = (X[17]-JVS[82]*X[69]-JVS[83]*X[70])/(JVS[81]);
  X[16] = (X[16]-JVS[79]*X[62]-JVS[80]*X[73])/(JVS[78]);
  X[15] = (X[15]-JVS[76]*X[64]-JVS[77]*X[69])/(JVS[75]);
  X[14] = (X[14]-JVS[73]*X[66]-JVS[74]*X[69])/(JVS[72]);
  X[13] = (X[13]-JVS[70]*X[69]-JVS[71]*X[72])/(JVS[69]);
  X[12] = (X[12]-JVS[67]*X[69]-JVS[68]*X[71])/(JVS[66]);
  X[11] = (X[11]-JVS[62]*X[23]-JVS[63]*X[48]-JVS[64]*X[61]-JVS[65]
         *X[73])/(JVS[61]);
  X[10] = (X[10]-JVS[60]*X[73])/(JVS[59]);
  X[9] = (X[9]-JVS[58]*X[61])/(JVS[57]);
  X[8] = (X[8]-JVS[56]*X[73])/(JVS[55]);
  X[7] = (X[7]-JVS[52]*X[27]-JVS[53]*X[37]-JVS[54]*X[69])/(JVS[51]);
  X[6] = (X[6]-JVS[49]*X[27]-JVS[50]*X[69])/(JVS[48]);
  X[5] = (X[5]-JVS[44]*X[62]-JVS[45]*X[64]-JVS[46]*X[66]-JVS[47]*X[72])
        /(JVS[43]);
  X[4] = (X[4]-JVS[41]*X[62]-JVS[42]*X[71])/(JVS[40]);
  X[3] = (X[3]-JVS[27]*X[46]-JVS[28]*X[48]-JVS[29]*X[50]-JVS[30]*X[51]
        -JVS[31]*X[52]-JVS[32]*X[61]-JVS[33]*X[62]-JVS[34]*X[63]
        -JVS[35]*X[64]-JVS[36]*X[65]-JVS[37]*X[66]-JVS[38]*X[67]
        -JVS[39]*X[72])/(JVS[26]);
  X[2] = (X[2]-JVS[18]*X[50]-JVS[19]*X[52]-JVS[20]*X[61]-JVS[21]*X[62]
        -JVS[22]*X[63]-JVS[23]*X[65]-JVS[24]*X[67]-JVS[25]*X[71])
        /(JVS[17]);
  X[1] = (X[1]-JVS[4]*X[19]-JVS[5]*X[26]-JVS[6]*X[43]-JVS[7]*X[46]
        -JVS[8]*X[48]-JVS[9]*X[49]-JVS[10]*X[50]-JVS[11]*X[51]-JVS[12]
        *X[52]-JVS[13]*X[53]-JVS[14]*X[61]-JVS[15]*X[68]-JVS[16]*X[73])
        /(JVS[3]);
  X[0] = (X[0]-JVS[1]*X[8]-JVS[2]*X[73])/(JVS[0]);
}

/* End of KppSolve function                                         */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* BLAS_UTIL - BLAS-LIKE utility functions                          */
/*   Arguments :                                                    */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/*--------------------------------------------------------------
 
  BLAS/LAPACK-like subroutines used by the integration algorithms
  It is recommended to replace them by calls to the optimized
       BLAS/LAPACK library for your machine
 
   (C) Adrian Sandu, Aug. 2004
 
--------------------------------------------------------------*/

#define ZERO  (double)0.0
#define ONE   (double)1.0
#define HALF  (double)0.5
#define TWO   (double)2.0
#define MOD(A,B) (int)((A)%(B))

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void WCOPY(int N, double X[], int incX, volatile double Y[], int incY)
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    copies a vector, x, to a vector, y:  y <- x
    only for incX=incY=1
    after BLAS
    replace this by the function from the optimized BLAS implementation:
        CALL  SCOPY(N,X,1,Y,1)   or   CALL  DCOPY(N,X,1,Y,1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
{
#if USE_SDK_BLAS == 1
    dcopy_spu(X, Y, N);
#else
      int i, M;
      if (N <= 0) return;

      M = MOD(N,8);
      if( M != 0 ) {
        for ( i = 0; i < M; i++ )
          Y[i] = X[i];
        if( N < 8 ) return;
      } 	
      for ( i = M; i<N; i+=8 ) {
        Y[i] = X[i];
        Y[i + 1] = X[i + 1];
        Y[i + 2] = X[i + 2];
        Y[i + 3] = X[i + 3];
        Y[i + 4] = X[i + 4];
        Y[i + 5] = X[i + 5];
        Y[i + 6] = X[i + 6];
        Y[i + 7] = X[i + 7];
      }
#endif
} /* end function WCOPY */


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void WAXPY(int N, double Alpha, double X[], int incX, double Y[], int incY )
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    constant times a vector plus a vector: y <- y + Alpha*x
    only for incX=incY=1
    after BLAS
    replace this by the function from the optimized BLAS implementation:
        CALL SAXPY(N,Alpha,X,1,Y,1) or  CALL DAXPY(N,Alpha,X,1,Y,1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
{
#if USE_SDK_BLAS == 1
    daxpy_spu(X, Y, Alpha, N);
#else
      int i, M;

      if (Alpha == ZERO) return;
      if (N  <=  0) return;

      M = MOD(N,4);
      if( M != 0 ) {
        for ( i = 0; i < M; i++ )
          Y[i] = Y[i] + Alpha*X[i];
        if ( N < 4 ) return;
      }
      
      for ( i = M; i < N; i += 4 ) {
        Y[i] = Y[i] + Alpha*X[i];
        Y[i + 1] = Y[i + 1] + Alpha*X[i + 1];
        Y[i + 2] = Y[i + 2] + Alpha*X[i + 2];
        Y[i + 3] = Y[i + 3] + Alpha*X[i + 3];
      }
#endif

} /* end function  WAXPY */



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void WSCAL(int N, double Alpha, double X[], int incX)
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    constant times a vector: x(1:N) <- Alpha*x(1:N) 
    only for incX=incY=1
    after BLAS
    replace this by the function from the optimized BLAS implementation:
        CALL SSCAL(N,Alpha,X,1) or  CALL DSCAL(N,Alpha,X,1)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
{
#if USE_SDK_BLAS == 1
    dscal_spu(X, Alpha, N);
#else
      int i, M;
      
      if (Alpha == ONE) return;
      if (N  <=  0) return;

      M = MOD(N,5);
      if( M  !=  0 ) {
        if (Alpha == (-ONE))
          for ( i = 0; i < M; i++ )  X[i] = -X[i];
        else {
	  if (Alpha == ZERO)
            for ( i = 0; i < M; i++ ) X[i] = ZERO;
          else
            for ( i = 0; i < M; i++ ) X[i] = Alpha*X[i];
        } // end else
        if( N < 5 ) return;
      } // end if
      
      if (Alpha == (-ONE))
        for ( i = M; i<N; i+=5 ) {
          X[i]     = -X[i];
          X[i + 1] = -X[i + 1];
          X[i + 2] = -X[i + 2];
          X[i + 3] = -X[i + 3];
          X[i + 4] = -X[i + 4];
        } // end for
      else {
        if (Alpha == ZERO)
          for ( i = M; i < N; i += 5 ) {
            X[i]     = ZERO;
            X[i + 1] = ZERO;
            X[i + 2] = ZERO;
            X[i + 3] = ZERO;
            X[i + 4] = ZERO;
          } // end for
        else
          for ( i = M; i < N; i += 5 ) {
            X[i]     = Alpha*X[i];
            X[i + 1] = Alpha*X[i + 1];
            X[i + 2] = Alpha*X[i + 2];
            X[i + 3] = Alpha*X[i + 3];
            X[i + 4] = Alpha*X[i + 4];
           } // end for
      }  // else
#endif
    
} /* end function WSCAL */

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
double WLAMCH( char C )
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    returns epsilon machine
    after LAPACK
    replace this by the function from the optimized LAPACK implementation:
         CALL SLAMCH('E') or CALL DLAMCH('E')
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
{
      int i;
      double Suma;
      static double Eps;
      static char First = 1;
      
      if (First) {
        First = 0;
        Eps = pow(HALF,16);
        for ( i = 17; i <= 80; i++ ) {
          Eps = Eps*HALF;
	  Suma = ONE + Eps;
	  if (Suma <= ONE) break;
        } /* end for */
        if (i==80) {
	   printf("\nERROR IN WLAMCH. Very small EPS = %g\n",Eps);
           return (double)2.2e-16;
	}
        Eps *= TWO; i--;      
      } /* end if First */

      return Eps;

} /* end function WLAMCH */


