/* Copyright (c) 2013, Intel Corporation
 *
 * Redistribution and use in source and binary forms, with or without 
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice, 
 *   this list of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, 
 *   this list of conditions and the following disclaimer in the documentation 
 *   and/or other materials provided with the distribution.
 * - Neither the name of Intel Corporation nor the names of its contributors 
 *   may be used to endorse or promote products derived from this software 
 *   without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */

// icpmain.cpp : entry point for the console application.
#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <sys/mman.h>
#include "dtime.h"

//
// Cache Block Tweaking Parameters
// 
unsigned int BlockA = 16;
unsigned int BlockB = 2560;
 
//
// Memory Routines
//
void *icpmalloc2(size_t numBytes)
{
void *ptr = NULL;
#if __MIC__
        ptr = mmap(0,numBytes,PROT_READ|PROT_WRITE,MAP_ANONYMOUS|MAP_PRIVATE|MAP_HUGETLB,-1,0);
#else
        ptr = mmap(0,numBytes,PROT_READ|PROT_WRITE,MAP_ANONYMOUS|MAP_PRIVATE,-1,0);
#endif
if( ptr == NULL ) printf("Problem calling mmap: ptr=NULL\r\n");
return(ptr);
}

void *icpmalloc(size_t numBytes)
{
size_t mic_alignment = 64;
void __attribute__((align(64))) *memptr = NULL;
int iret = posix_memalign( &memptr, mic_alignment, numBytes );
if( iret != 0 ) printf("problem calling posix_memalign: %d\r\n",iret);
return( memptr );
}

void icpfree2(void *ptr)
{
	// munmap( ptr, len ); // don't have len available normally
}
void icpfree(void *ptr)
{
	free(ptr);
}

#if defined(AOS)
#define AOSFLAG     1              // 1 for AOS
#else
#define AOSFLAG     0              // 0 for SOA
#endif

#if defined(DOUBLEPREC)
#define FPPRECFLAG  2              // 2 for Double precision
#define FLOATINGPTPRECISION double
#else
#define FPPRECFLAG  1              // 1 for Single precision
#define FLOATINGPTPRECISION float
#endif

#define DEBUGOUTPUT 0              // 0 for no output, 1 for debug output
//#define UNSIGNEDTYPE  unsigned 
#define UNSIGNEDTYPE    

// 
// Structure to do Array of Structures tests (vs. Structures or Set of Arrays)
typedef struct Point3d 
{
	//FLOATINGPTPRECISION x,y,z,t; // padding for 4x multiples 
	FLOATINGPTPRECISION x,y,z;    // no padding option
}  Point3d,*Point3dPtr;

#define RETAIN free_if(0)
#define REUSE  alloc_if(0)
#define ALLOC  alloc_if(1)
#define FREE   free_if(1)

int  eigend(double aptr[],double rptr[],int nsize);
void qcvcrs(double qmatrx[][4],FLOATINGPTPRECISION covcrs[][3]);
void qat2rot(FLOATINGPTPRECISION rot[3][3], FLOATINGPTPRECISION q[4]);
void listAccumulatedValues(FILE *fptext,
			   FLOATINGPTPRECISION ttotresid,	
			   FLOATINGPTPRECISION tavgxfm[3],
			   FLOATINGPTPRECISION tavgorg[3],
			   FLOATINGPTPRECISION tcrosscov[3][3]);
#if AOSFLAG == 1 
void listDataArrays(FILE *fptext, Point3dPtr org, Point3dPtr tfm,
		    int nptsmod16, int num2show);
#endif
#if AOSFLAG == 0 
void listDataArrays(FILE *fptext, FLOATINGPTPRECISION *orgx, FLOATINGPTPRECISION *tfmx,
		    int nptsmod16, int num2show);
#endif

void qcvcrs(double qmatrx[4][4],FLOATINGPTPRECISION covcrs[3][3])
{
/*--------------------------------------------------------------------*/
	qmatrx[0][0] = covcrs[0][0] + covcrs[1][1] + covcrs[2][2];
	qmatrx[0][1] = covcrs[1][2] - covcrs[2][1];
	qmatrx[1][0] = qmatrx[0][1];
	qmatrx[0][2] = covcrs[2][0] - covcrs[0][2];
	qmatrx[2][0] = qmatrx[0][2];
	qmatrx[0][3] = covcrs[0][1] - covcrs[1][0];
	qmatrx[3][0] = qmatrx[0][3];
	qmatrx[1][1] = covcrs[0][0] + covcrs[0][0] - qmatrx[0][0];
	qmatrx[2][2] = covcrs[1][1] + covcrs[1][1] - qmatrx[0][0];
	qmatrx[3][3] = covcrs[2][2] + covcrs[2][2] - qmatrx[0][0];
	qmatrx[2][1] = covcrs[0][1] + covcrs[1][0];
	qmatrx[1][2] = qmatrx[2][1];
	qmatrx[3][2] = covcrs[1][2] + covcrs[2][1];
	qmatrx[2][3] = qmatrx[3][2];
	qmatrx[1][3] = covcrs[2][0] + covcrs[0][2];
	qmatrx[3][1] = qmatrx[1][3];
	return;
/*--------------------------------------------------------------------*/
}

void qat2rot(FLOATINGPTPRECISION rot[3][3], FLOATINGPTPRECISION q[4])
{
	double q0,q1,q2,q3;
/*--------------------------------------------------------------------*/
	q0 = (double)q[0]; q1 = (double)q[1]; 
 	q2 = (double)q[2]; q3 = (double)q[3];
	rot[0][0] = (FLOATINGPTPRECISION)(q0*q0 + q1*q1 - ( q2*q2 + q3*q3 ));
	rot[0][1] = (FLOATINGPTPRECISION)(2.0*( q1*q2 - q0*q3 ));
	rot[0][2] = (FLOATINGPTPRECISION)(2.0*( q1*q3 + q0*q2 ));
	rot[1][0] = (FLOATINGPTPRECISION)(2.0*( q1*q2 + q0*q3 ));
	rot[1][1] = (FLOATINGPTPRECISION)(q0*q0 + q2*q2 - ( q3*q3 + q1*q1 ));
	rot[1][2] = (FLOATINGPTPRECISION)(2.0*( q2*q3 - q0*q1 ));
	rot[2][0] = (FLOATINGPTPRECISION)(2.0*( q1*q3 - q0*q2 ));
	rot[2][1] = (FLOATINGPTPRECISION)(2.0*( q2*q3 + q0*q1 ));
	rot[2][2] = (FLOATINGPTPRECISION)(q0*q0 + q3*q3 - ( q1*q1 + q2*q2 ));
	return;
/*--------------------------------------------------------------------*/
}

#define MAXTHREADS 256

int main(int argc, char* argv[])
{
	UNSIGNEDTYPE int ithread = 0;
	UNSIGNEDTYPE int max_hostnative_threads = 1;
	UNSIGNEDTYPE int num_mic_threads = 0;
	UNSIGNEDTYPE int num_host_threads = 0;
	UNSIGNEDTYPE int nargs = 0;
	UNSIGNEDTYPE int nu,nv;
	UNSIGNEDTYPE int iu,iv;
	UNSIGNEDTYPE int i,j;
	UNSIGNEDTYPE int norgpts,nxfmpts;
	UNSIGNEDTYPE int micBlockA,micBlockB;
	UNSIGNEDTYPE int niter;
	UNSIGNEDTYPE int select;
	UNSIGNEDTYPE int nptsmod16 = 0;

#if AOSFLAG == 1
	Point3dPtr org = NULL;
	Point3dPtr tfm = NULL;
#endif 
#if AOSFLAG == 0 
	FLOATINGPTPRECISION __attribute__((align(64))) *orgx = NULL;
	FLOATINGPTPRECISION __attribute__((align(64))) *orgy = NULL;
	FLOATINGPTPRECISION __attribute__((align(64))) *orgz = NULL;
	FLOATINGPTPRECISION __attribute__((align(64))) *tfmx = NULL;
	FLOATINGPTPRECISION __attribute__((align(64))) *tfmy = NULL;
	FLOATINGPTPRECISION __attribute__((align(64))) *tfmz = NULL;
#endif
        FLOATINGPTPRECISION q0,q1,q2,q3;
	FLOATINGPTPRECISION tx,ty,tz;
	FLOATINGPTPRECISION x,y,z;
	FLOATINGPTPRECISION ru,rv;
	FLOATINGPTPRECISION Rinit[3][3],Tinit[3],qinit[4];
	FLOATINGPTPRECISION Rf[3][3],Tf[3],qf[4];
	FLOATINGPTPRECISION qnorm;
	FLOATINGPTPRECISION factor;
	FLOATINGPTPRECISION packedvals[MAXTHREADS][2][4][4];
	FLOATINGPTPRECISION tcrosscov[3][3];
	FLOATINGPTPRECISION tavgxfm[3];
	FLOATINGPTPRECISION tavgorg[3];
	FLOATINGPTPRECISION firstresid;
	FLOATINGPTPRECISION ttotresid;
	FLOATINGPTPRECISION lasttotresid;
	FLOATINGPTPRECISION valueOfOne = 1.0;
	FLOATINGPTPRECISION valueOfHalf = 0.5;
	FLOATINGPTPRECISION valueOfTwo = 2.0;
	FLOATINGPTPRECISION valueOfZero = 0.0;
	FLOATINGPTPRECISION *minkDist2 = NULL;

	// always doubles no matter where points are single or double
	double eveca[4][4];
	double emat[4][4];
	double emax;
	double totaltime = valueOfZero;
	double titer = valueOfZero;
	double totalestimate = valueOfZero;
	double gigaflops = valueOfZero;
	double tgigaflops = valueOfZero;
	double tlooptime = valueOfZero;

	unsigned int mod16,mod32,mod64;
	int num2show = 10;
	int hostOffloadNative = 0; // 0=host // 1=offload // 2=native
	int uiNumThreads = 0;
	int *minkLoc = NULL;
//==========================================================//
 	dtimeInit();
	totaltime = dtimeGet();
	nu = nv = 500; q0 = valueOfOne; 
	tx = 0.3;
	ty = 0.4;
	tz = 0.5;
	q1 = q2 = q3 = 0.1;
	if( argc < 5 )
	{
		printf("Usage: %s {0=host|1=off|2=native} nthreads nu nv [q0 q1 q2 q3 tx ty tz]\r\n",argv[0]);
		printf("Notes: numThreads=24/32 for Host, N-4 for Offload, N for Native\r\n");
		printf("Notes: Runtime increases as 4th power of nu=nv parameter.\r\n");
		exit(-1);
	}
	printf("Iterative Closest Point (ICP) Benchmark Release 10.0 (for B0/B1/Product KNC) (Cache-Blocked)\r\n");
	printf("Single Source for Host(OpenMP,Serial), MIC-Native(OpenMP,Serial) & MIC-Offload Executables\r\n");
	//
	// parse command line values
	// 
	if( argc >  1 ) nargs = sscanf(argv[1],"%d",&hostOffloadNative); 
	if( strstr(argv[0],"host") )    hostOffloadNative = 0;
	if( strstr(argv[0],"offload") ) hostOffloadNative = 1;
	if( strstr(argv[0],"mic") )     hostOffloadNative = 2;
	if( argc >  2 ) nargs = sscanf(argv[2],"%d",&uiNumThreads);
	if( argc >  3 ) nargs = sscanf(argv[3],"%d",&nu); 
	if( argc >  4 ) nargs = sscanf(argv[4],"%d",&nv);
#if FPPRECFLAG == 1
	if( argc >  5 ) nargs = sscanf(argv[5],"%f",&q0); 
	if( argc >  6 ) nargs = sscanf(argv[6],"%f",&q1);
	if( argc >  7 ) nargs = sscanf(argv[7],"%f",&q2); 
	if( argc >  8 ) nargs = sscanf(argv[8],"%f",&q3);
	if( argc >  9 ) nargs = sscanf(argv[9],"%f",&tx); 
	if( argc > 10 ) nargs = sscanf(argv[10],"%f",&ty);
	if( argc > 11 ) nargs = sscanf(argv[11],"%f",&tz);
#elif FPPRECFLAG == 2
	if( argc >  5 ) nargs = sscanf(argv[5],"%lf",&q0); 
	if( argc >  6 ) nargs = sscanf(argv[6],"%lf",&q1);
	if( argc >  7 ) nargs = sscanf(argv[7],"%lf",&q2); 
	if( argc >  8 ) nargs = sscanf(argv[8],"%lf",&q3);
	if( argc >  9 ) nargs = sscanf(argv[9],"%lf",&tx); 
	if( argc > 10 ) nargs = sscanf(argv[10],"%lf",&ty);
	if( argc > 11 ) nargs = sscanf(argv[11],"%lf",&tz);
#endif
	printf("Host|Offload|Native: %d  (0=Host,1=Offload,2=Native)\r\n",hostOffloadNative);
	printf("Number_Of_Threads: %d\r\n",uiNumThreads);
	printf("q: %g %g %g %g / t: %g %g %g\r\n",q0,q1,q2,q3,tx,ty,tz);
	printf("nu,nv= %d %d\r\n",nu,nv);
	//
	// make a rotation matrix after normalizing the quaternion
	// save initial state of everything
	//
	if( q0 < 0.0 ) { q0 = -q0; q1 = -q1; q2 = -q2; q3 = -q3; }
	qnorm = q0*q0 + q1*q1 + q2*q2 + q3*q3;
	qnorm = valueOfOne/sqrt( qnorm );
	q0 *= qnorm; q1 *= qnorm; q2 *= qnorm; q3 *= qnorm;
	qinit[0] = q0; qinit[1] = q1; qinit[2] = q2; qinit[3] = q3;
	qat2rot(Rinit,qinit);
	memcpy(Rf,Rinit,sizeof(FLOATINGPTPRECISION)*9);
	Tinit[0] = tx; Tinit[1] = ty; Tinit[2] = tz;
	memcpy(Tf,Tinit,sizeof(FLOATINGPTPRECISION)*3);
	memcpy(qf,qinit,sizeof(FLOATINGPTPRECISION)*4);
	//
	// Create 2 point sets: one original (org) and one transformed (tfm)
	//
	micBlockA = BlockA;
	micBlockB = BlockB;
	nxfmpts = norgpts = nu*nv;
	nptsmod16 = (1+(norgpts>>4))<<4;

	minkLoc = (int *)malloc(sizeof(int)*nptsmod16);
	minkDist2 = (FLOATINGPTPRECISION *)malloc(sizeof(FLOATINGPTPRECISION)*nptsmod16);
	memset(minkLoc,-1,sizeof(int)*nptsmod16);
	for(i=0;i<nptsmod16;i++) minkDist2[i] = 1.0e38f;

#if AOSFLAG == 1
	org = (Point3dPtr)icpmalloc(sizeof(Point3d)*2*nptsmod16);
	tfm = org + nptsmod16;
	if( org == NULL ) 
	{
		printf("Point array Memory allocation on Host failed!\r\n"); exit(-1); 
	}
	memset(org,0,sizeof(Point3d)*2*nptsmod16);
        //
        // check & report alignment problem
        mod16 = ((unsigned long long)org%16);
        mod32 = ((unsigned long long)org%32);
        mod64 = ((unsigned long long)org%64);
	printf("mod64,mod32,mod16= %d %d %d (should all be zero)\r\n",mod64,mod32,mod16);
#endif
#if AOSFLAG == 0
	orgx = (FLOATINGPTPRECISION *)icpmalloc(sizeof(FLOATINGPTPRECISION)*6*nptsmod16);
	orgy = orgx + nptsmod16; orgz = orgy + nptsmod16;
	tfmx = orgx + 3*nptsmod16;
	tfmy = tfmx + nptsmod16; tfmz = tfmy + nptsmod16;
	if( orgx == NULL ) 
	{
		printf("Point array Memory allocation on Host failed!\r\n"); exit(-1); 
	}
	memset(orgx,0,sizeof(FLOATINGPTPRECISION)*6*nptsmod16);
        //
        // check alignment alignment problem
        mod16 = ((unsigned long long)orgx%16);
        mod32 = ((unsigned long long)orgx%32);
        mod64 = ((unsigned long long)orgx%64);
	printf("mod64,mod32,mod16= %d %d %d (should all be zero)\r\n",mod64,mod32,mod16);
#endif
	ru = valueOfOne/(FLOATINGPTPRECISION)nu; rv = valueOfOne/(FLOATINGPTPRECISION)nv;
	//
	// Create Original and Xformed Data Sets 
	//
	for(iv=0;iv<nv;iv++)
	{
		y = (FLOATINGPTPRECISION)iv;
		y = valueOfTwo*(ru*y - valueOfHalf);

		for(iu=0;iu<nu;iu++)
		{
			x = (FLOATINGPTPRECISION)iu;
			x = valueOfTwo*(ru*x - valueOfHalf);
			z = valueOfHalf*(x*x + valueOfHalf*y*y);
			i = nu*iv + iu;

#if AOSFLAG == 1
			org[i].x = x; org[i].y = y; org[i].z = z;
			tfm[i].x = Rf[0][0]*x + Rf[0][1]*y + Rf[0][2]*z + Tf[0];
			tfm[i].y = Rf[1][0]*x + Rf[1][1]*y + Rf[1][2]*z + Tf[1];
			tfm[i].z = Rf[2][0]*x + Rf[2][1]*y + Rf[2][2]*z + Tf[2];
#endif
#if AOSFLAG == 0
			orgx[i] = x; orgy[i] = y; orgz[i] = z;
			tfmx[i] = Rf[0][0]*x + Rf[0][1]*y + Rf[0][2]*z + Tf[0];
			tfmy[i] = Rf[1][0]*x + Rf[1][1]*y + Rf[1][2]*z + Tf[1];
			tfmz[i] = Rf[2][0]*x + Rf[2][1]*y + Rf[2][2]*z + Tf[2];
#endif
		}
	}
		
	num_mic_threads = 0;
#ifdef _OPENMP
	max_hostnative_threads = omp_get_max_threads();
#endif 
	printf("Max {Host|Native} Threads: %d\r\n",max_hostnative_threads);
	if( hostOffloadNative == 1 )  // Offload Case
	{
		num_mic_threads = uiNumThreads; 
#ifdef _OPENMP
		num_host_threads = 2;
#pragma omp parallel
{
#pragma omp single 
{
	        omp_set_num_threads(num_host_threads); // use 2 threads on host
}
}
#endif 
	}
	else if( (hostOffloadNative&(0x1)) == 0 )  // Even Value Native or Host Case
	{
		num_host_threads = num_mic_threads = uiNumThreads;
#ifdef _OPENMP
	        omp_set_num_threads(num_mic_threads); // on mic native or host
#endif 
	}
//
// copy xyzorg over one time to the card
// 
#if AOSFLAG == 1
#pragma offload target(mic:0) in(org:length(nptsmod16) RETAIN) 
{ ; }
#endif
#if AOSFLAG == 0
#pragma offload target(mic:0) in(orgx:length(3*nptsmod16) RETAIN) 
{ ; }
#endif
//
// set number of threads: should be in parallel region to change
// should be in single region in parallel region to avoid thrashing all the cores 
//
#ifdef _OPENMP
#pragma offload target(mic:0) inout(num_mic_threads,nargs)
{
        printf(" OpenMP: Serial Region: ThreadId= %2d | omp_max_threads= %2d | omp_num_threads= %2d\r\n",
                omp_get_thread_num(), omp_get_max_threads(), omp_get_num_threads());

#pragma omp parallel 
        {
#pragma omp single
            {
	        omp_set_num_threads(num_mic_threads); // done in single part of parallel region 
                printf(" OpenMP: Single Region in Parallel Region: ThreadId= %2d | omp_max_threads= %2d | omp_num_threads= %2d\r\n",
                        omp_get_thread_num(), omp_get_max_threads(),omp_get_num_threads() );
            }
         }
} // end of pragma offload region
#endif

#if DEBUGOUTPUT == 1

#if AOSFLAG == 1
        listDataArrays(stdout,org,tfm,nptsmod16,num2show);
#endif
#if AOSFLAG == 0
        listDataArrays(stdout,orgx,tfmx,nptsmod16,num2show);
#endif

#endif


	printf("NumComputeThreads: %d\r\n",num_mic_threads);
	if( hostOffloadNative == 1 ) printf("NumHostThreads: %d\r\n",num_host_threads);
	fflush(stdout);
	//
        // Bring two data sets back together using ICP algorithm
	//
	niter = 0; lasttotresid = 1.0e38f;
	tlooptime = dtimeGet();
	while(1)
	{
		titer = dtimeGet();
		// for each point on the xform'd object find closest point on original object
		// accumulate cross covariance matrix and average vectors (centroids)
		//
		memset(tcrosscov,0,sizeof(FLOATINGPTPRECISION)*9); ttotresid = valueOfZero;
		memset(tavgxfm,0,sizeof(FLOATINGPTPRECISION)*3);   memset(tavgorg,0,sizeof(FLOATINGPTPRECISION)*3);
		memset(packedvals,0,2*sizeof(FLOATINGPTPRECISION)*16*MAXTHREADS);
	        memset(minkLoc,-1,sizeof(int)*nptsmod16);
		for(i=0;i<nptsmod16;i++) minkDist2[i] = 1.0e38f;
		//
		// Do MIC offload of primary compute loop
		//
#if AOSFLAG == 1
#pragma offload target(mic:0) nocopy(org:length(nptsmod16) REUSE RETAIN) \
in(tfm:length(nptsmod16) ) \
in(minkLoc:length(nptsmod16) ) \
in(minkDist2:length(nptsmod16) ) \
inout(packedvals) \
in(nxfmpts,norgpts,nptsmod16,micBlockA,micBlockB)
#endif
#if AOSFLAG == 0
#pragma offload target(mic:0) nocopy(orgx:length(3*nptsmod16) REUSE RETAIN) \
in(tfmx:length(3*nptsmod16) ) \
in(minkLoc:length(nptsmod16) ) \
in(minkDist2:length(nptsmod16) ) \
inout(packedvals) \
in(nxfmpts,norgpts,nptsmod16,micBlockA,micBlockB)
#endif
{
//==============================================================================
#if AOSFLAG == 1
#pragma omp parallel for shared(org,tfm,packedvals,nxfmpts,norgpts) 
#endif
#if AOSFLAG == 0
#pragma omp parallel for shared(orgx,tfmx,packedvals,nxfmpts,norgpts) 
#endif
		for(unsigned int kblock=0;kblock<nxfmpts;kblock+=micBlockA)
		{
		   unsigned int imicthread;
#if AOSFLAG == 0
		   FLOATINGPTPRECISION *micorgy;
		   FLOATINGPTPRECISION *micorgz;
		   FLOATINGPTPRECISION *mictfmy;
		   FLOATINGPTPRECISION *mictfmz;
		   micorgy = orgx + nptsmod16; micorgz = micorgy + nptsmod16;
		   mictfmy = tfmx + nptsmod16; mictfmz = mictfmy + nptsmod16;
#endif
#ifdef _OPENMP
		   imicthread = omp_get_thread_num();
#endif
#ifndef _OPENMP
		   imicthread = 0;
#endif
		   unsigned int kmax = kblock+micBlockA;
		   if( kmax > nxfmpts ) kmax = nxfmpts;
 		   for(unsigned int mblock=0;mblock<norgpts;mblock+=micBlockB)
		   {
		       unsigned int mmax = mblock+micBlockB;
		       if( mmax > norgpts ) mmax = norgpts;
		       for(unsigned int k=kblock;k<kmax;k++)
	  	       {
			   FLOATINGPTPRECISION micdist2;
			   FLOATINGPTPRECISION micmindist2;
			   FLOATINGPTPRECISION xmic,ymic,zmic;
   			   int micminjloc;

#if AOSFLAG == 1
			   xmic = tfm[k].x; ymic = tfm[k].y; zmic = tfm[k].z;
#endif
#if AOSFLAG == 0
			   xmic = tfmx[k]; ymic = mictfmy[k]; zmic = mictfmz[k];
#endif

			   micminjloc = minkLoc[k]; 
			   micmindist2 = minkDist2[k];

#pragma vector aligned
#ifdef __MIC__
#pragma unroll(0)
#endif
#ifndef __MIC__
#pragma unroll(4)
#endif
			   for(unsigned int m=mblock;m<mmax;++m)
			   {
			        FLOATINGPTPRECISION micdx,micdy,micdz;
#if AOSFLAG == 1
				micdx =    tfm[k].x -    org[m].x;
				micdy =    tfm[k].y -    org[m].y;
				micdz =    tfm[k].z -    org[m].z;
#endif
#if AOSFLAG == 0
				micdx =    xmic -    orgx[m];
				micdy =    ymic - micorgy[m];
				micdz =    zmic - micorgz[m];
#endif
				micdist2 = micdx*micdx + micdy*micdy + micdz*micdz;
				if( micdist2 < micmindist2 )
				{
					micmindist2 = micdist2;
					micminjloc = m;
				}
			   } //m
			   minkLoc[k] = micminjloc;
			   minkDist2[k] = micmindist2;
			} //k
	           } //mblock 


		   for(unsigned int k=kblock;k<kmax;k++)
	  	   {
			FLOATINGPTPRECISION xmic,ymic,zmic;
			FLOATINGPTPRECISION micxorg, micyorg, miczorg;
#if AOSFLAG == 1
			xmic = tfm[k].x; ymic = tfm[k].y; zmic = tfm[k].z;
#endif
#if AOSFLAG == 0
			xmic = tfmx[k]; ymic = mictfmy[k]; zmic = mictfmz[k];
#endif
			packedvals[imicthread][0][0][3] += xmic;
			packedvals[imicthread][0][1][3] += ymic;
			packedvals[imicthread][0][2][3] += zmic;
			packedvals[imicthread][0][3][3] += minkDist2[k];
			if( minkLoc[k] >= 0 )
			{
#if AOSFLAG == 1
				micxorg = org[minkLoc[k]].x;
				micyorg = org[minkLoc[k]].y;
				miczorg = org[minkLoc[k]].z;
#endif 
#if AOSFLAG == 0
				micxorg = orgx[minkLoc[k]];
				micyorg = micorgy[minkLoc[k]];
				miczorg = micorgz[minkLoc[k]];
#endif 
				packedvals[imicthread][0][3][0] += micxorg;
				packedvals[imicthread][0][3][1] += micyorg;
				packedvals[imicthread][0][3][2] += miczorg;

				packedvals[imicthread][0][0][0] += (micxorg*xmic);
				packedvals[imicthread][0][0][1] += (micyorg*xmic);
				packedvals[imicthread][0][0][2] += (miczorg*xmic);
				packedvals[imicthread][0][1][0] += (micxorg*ymic);
				packedvals[imicthread][0][1][1] += (micyorg*ymic);
				packedvals[imicthread][0][1][2] += (miczorg*ymic);
				packedvals[imicthread][0][2][0] += (micxorg*zmic);
				packedvals[imicthread][0][2][1] += (micyorg*zmic);
				packedvals[imicthread][0][2][2] += (miczorg*zmic);
			}
	  	   } //k
		}// kblock // an ithread is assigned to each
//==============================================================================
}// end of mic offload section
		//
		// do explicit accumulation/reduction on host
		//
	        for(ithread=0;ithread<num_mic_threads;++ithread)
 		{
			tavgxfm[0] += packedvals[ithread][0][0][3];
			tavgxfm[1] += packedvals[ithread][0][1][3];
			tavgxfm[2] += packedvals[ithread][0][2][3];
			tavgorg[0] += packedvals[ithread][0][3][0];
			tavgorg[1] += packedvals[ithread][0][3][1];
			tavgorg[2] += packedvals[ithread][0][3][2];
			tcrosscov[0][0] += packedvals[ithread][0][0][0];
			tcrosscov[0][1] += packedvals[ithread][0][0][1];
			tcrosscov[0][2] += packedvals[ithread][0][0][2];
			tcrosscov[1][0] += packedvals[ithread][0][1][0];
			tcrosscov[1][1] += packedvals[ithread][0][1][1];
			tcrosscov[1][2] += packedvals[ithread][0][1][2];
			tcrosscov[2][0] += packedvals[ithread][0][2][0];
			tcrosscov[2][1] += packedvals[ithread][0][2][1];
			tcrosscov[2][2] += packedvals[ithread][0][2][2];
			ttotresid += packedvals[ithread][0][3][3];
		}
//==============================================================================
	        titer = dtimeGet() - titer;
                //
                //  back on host with small 3x3 data in offload version
                // 
		factor = valueOfOne/(FLOATINGPTPRECISION)nxfmpts;
		ttotresid *= (factor);
		for(i=0;i<3;i++)
		{
			tavgxfm[i] *= (factor);
			tavgorg[i] *= (factor);
			for(j=0;j<3;j++) { tcrosscov[i][j] *= (factor); }
		}
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				tcrosscov[i][j] -= (tavgxfm[i]*tavgorg[j]);
			}
		}
		//
		//	optionally list values to stdout
		//
#if DEBUGOUTPUT == 1
		listAccumulatedValues(stdout,ttotresid,tavgxfm,tavgorg,tcrosscov);
#endif
		//
		// convert from 3x3 unsymmetric cross covariance to 4x4 symmetric 
		//
		qcvcrs(emat,tcrosscov); 
		eigend((double *)emat,(double *)eveca,4);

		select = 0;
		emax = emat[0][0];
		if( emax < emat[1][1] ) { emax = emat[1][1]; select = 1; }
		if( emax < emat[2][2] ) { emax = emat[2][2]; select = 2; }
		if( emax < emat[3][3] ) { emax = emat[3][3]; select = 3; }
		q0 = (FLOATINGPTPRECISION)eveca[0][select]; 
		q1 = (FLOATINGPTPRECISION)eveca[1][select];
		q2 = (FLOATINGPTPRECISION)eveca[2][select]; 
		q3 = (FLOATINGPTPRECISION)eveca[3][select];
		if( q0 < valueOfZero ) { q0 = -q0; q1 = -q1; q2 = -q2; q3 = -q3; }
		qnorm = q0*q0 + q1*q1 + q2*q2 + q3*q3;
		qnorm = valueOfOne/sqrt( qnorm );
		q0 *= qnorm; q1 *= qnorm; q2 *= qnorm; q3 *= qnorm;
		qf[0] = q0; qf[1] = q1; qf[2] = q2; qf[3] = q3;
		qat2rot(Rf,qf);

		Tf[0] = Rf[0][0]*tavgxfm[0] + Rf[0][1]*tavgxfm[1] + Rf[0][2]*tavgxfm[2];
		Tf[1] = Rf[1][0]*tavgxfm[0] + Rf[1][1]*tavgxfm[1] + Rf[1][2]*tavgxfm[2];
		Tf[2] = Rf[2][0]*tavgxfm[0] + Rf[2][1]*tavgxfm[1] + Rf[2][2]*tavgxfm[2];
		Tf[0] = tavgorg[0] - Tf[0];
		Tf[1] = tavgorg[1] - Tf[1];
		Tf[2] = tavgorg[2] - Tf[2];
		//
		// Update xfm point locations for now (not best way to do it numerically)
		//
//==============================================================================
#if AOSFLAG == 1
#pragma omp parallel for shared(tfm,Rf,Tf,nxfmpts) private(x,y,z)
		for(i=0;i<nxfmpts;i++)
		{
			x = tfm[i].x; y = tfm[i].y; z = tfm[i].z;
			tfm[i].x = Rf[0][0]*x + Rf[0][1]*y + Rf[0][2]*z + Tf[0];
			tfm[i].y = Rf[1][0]*x + Rf[1][1]*y + Rf[1][2]*z + Tf[1];
			tfm[i].z = Rf[2][0]*x + Rf[2][1]*y + Rf[2][2]*z + Tf[2];
		}
#endif
#if AOSFLAG == 0
#pragma omp parallel for shared(tfmx,tfmy,tfmz,Rf,Tf,nxfmpts) private(x,y,z)
		for(i=0;i<nxfmpts;i++)
		{
			x = tfmx[i]; y = tfmy[i]; z = tfmz[i];
			tfmx[i] = Rf[0][0]*x + Rf[0][1]*y + Rf[0][2]*z + Tf[0];
			tfmy[i] = Rf[1][0]*x + Rf[1][1]*y + Rf[1][2]*z + Tf[1];
			tfmz[i] = Rf[2][0]*x + Rf[2][1]*y + Rf[2][2]*z + Tf[2];
		}
#endif
//==============================================================================

		if( fabs(lasttotresid - ttotresid) < 1.0e-4f ) 
		{
		        printf("%d: ttotresid: %g\r\n",niter,ttotresid);
		        printf("%d: lasttotresid: %g\r\n",niter,lasttotresid);
		        printf("%d: diff: %g\r\n",niter, fabs(lasttotresid - ttotresid) );
			break;
		}
		printf("%d: ttotresid: %g (dt= %.3f)\r\n",niter,ttotresid,titer);
		fflush(stdout);
		if( niter <= 0 ) firstresid = ttotresid;
		lasttotresid = ttotresid;
		++niter;
		tgigaflops += titer;
	}
	tlooptime = dtimeGet() - tlooptime;
	gigaflops = 8.0*(double)nu*nu*nv*nv*niter;
	if( tgigaflops > 0.0 ) gigaflops /= tgigaflops;

	printf("Total Iterations: %d\r\n",niter);
	printf("<<< GigaFlops: %.3f >>>\r\n",1.0e-9*gigaflops);
	totaltime = dtimeSince(totaltime,"Total ICP Runtime");
	{
		long tsecs = time(NULL);
		//printf("tsecs= %d\r\n",tsecs);
		char filename[4096];
		sprintf(filename,"icpruntime%d.log",(int)tsecs-1333000000);
		FILE *fp = fopen(filename,"wb");
		if( fp )
		{
		    fprintf(fp,"Iterative Closest Point Benchmark v6\r\n");
		    fprintf(fp,"HostOffloadNative= %d (0=host,1=offload,2=native\r\n",hostOffloadNative);
		    fprintf(fp,"Runtime= %g seconds\r\n",totaltime);
		    fprintf(fp,">>>GigaFlops: %.3f\r\n",1.0e-9*gigaflops);
		    fprintf(fp,"%d: ttotresid: %g\r\n",0,firstresid);
		    fprintf(fp,"%d: ttotresid: %g\r\n",niter,ttotresid);
		    for(i=0;i<argc;i++)
		    { fprintf(fp,"argv[%d]= '%s'\r\n",i,argv[i]);
		    }
		    fprintf(fp,"%d: ttotresid: %g\r\n",niter,ttotresid);
		    fprintf(fp,"Runtime= %g seconds\r\n",totaltime);
#if AOSFLAG == 1
		    fprintf(fp,"Array of Structures (AoS) Compilation\r\n");
#endif
#if AOSFLAG == 0
		    fprintf(fp,"Set of Arrays (SoA) Compilation (aka Structure of Arrays)\r\n");
#endif
#if FPPRECFLAG == 1
		    fprintf(fp,"Single Precision (SP) Compilation\r\n");
#endif
#if FPPRECFLAG == 2
		    fprintf(fp,"Double Precision (DP) Compilation\r\n");
#endif
	            //if( getenv("MIC_OMP_NUM_THREADS") != NULL ) 
		    //    fprintf(fp,"MIC_OMP_NUM_THREADS= %s'\r\n", getenv("MIC_OMP_NUM_THREADS"));
	            //if( getenv("OMP_NUM_THREADS") != NULL ) 
		    //    fprintf(fp,"OMP_NUM_THREADS= %s'\r\n", getenv("OMP_NUM_THREADS"));
		    fprintf(fp,"NumThreads= %d\r\n",uiNumThreads);
		    fprintf(fp,"MicThreads= %d\r\n",num_mic_threads);
		    fclose(fp); fp = NULL;
		}
		sprintf(filename,"icp.csv");
		fp = fopen(filename,"ab+");
		if( fp )
		{
			fprintf(fp,"icpv8, %d, %d, %s, %d, %d, %d, %d, %g, %.3f\r\n",
				hostOffloadNative, 
				(int)tsecs-1333000000, // time stamp
				argv[0],nu,nv,uiNumThreads,num_mic_threads,
				totaltime,1.0e-9*gigaflops);
		        fclose(fp); fp = NULL;
		}
	}
#if AOSFLAG == 1
	icpfree(org);
#endif
#if AOSFLAG == 0
	icpfree(orgx);
#endif
	return 0;
}


#if AOSFLAG == 1 
void listDataArrays(FILE *fptext,
		    Point3dPtr org,
		    Point3dPtr tfm,
		    int nptsmod16, int num2show)
{
	fprintf(fptext,"NumPts= %d / Num2Show= %d\r\n",nptsmod16,num2show);
	for(int i=0;i<num2show;i++)
	{
		fprintf(fptext,"OrgXyz= %10.3f %10.3f %10.3f / TfmXyz= %10.3f  %10.3f  %10.3f\r\n",
			org[i].x,org[i].y,org[i].z,tfm[i].x,tfm[i].y,tfm[i].z);
	}
	return;
}
#endif

#if AOSFLAG == 0 
void listDataArrays(FILE *fptext,
		    FLOATINGPTPRECISION *orgx,
		    FLOATINGPTPRECISION *tfmx,
		    int nptsmod16, int num2show)
{
	FLOATINGPTPRECISION *orgy;
	FLOATINGPTPRECISION *orgz;
	FLOATINGPTPRECISION *tfmy;
	FLOATINGPTPRECISION *tfmz;
	orgy = orgx + nptsmod16; orgz = orgy + nptsmod16;
	tfmy = tfmx + nptsmod16; tfmz = tfmy + nptsmod16;
	fprintf(fptext,"NumPts= %d / Num2Show= %d\r\n",nptsmod16,num2show);
	for(int i=0;i<num2show;i++)
	{
		fprintf(fptext,"OrgXyz= %10.3f %10.3f %10.3f / TfmXyz= %10.3f  %10.3f  %10.3f\r\n",
			orgx[i],orgy[i],orgz[i],tfmx[i],tfmy[i],tfmz[i]);
	}
	return;
}
#endif


void listAccumulatedValues(FILE *fptext,
			   FLOATINGPTPRECISION ttotresid,	
			   FLOATINGPTPRECISION tavgxfm[3],
			   FLOATINGPTPRECISION tavgorg[3],
			   FLOATINGPTPRECISION tcrosscov[3][3])
{
	fprintf(fptext,"ttotresid: %g\r\n",ttotresid);
	fprintf(fptext,"tavgxfm: %g %g %g\r\n",tavgxfm[0],tavgxfm[1],tavgxfm[2]);
	fprintf(fptext,"tavgorg: %g %g %g\r\n",tavgorg[0],tavgorg[1],tavgorg[2]);
	fprintf(fptext,"tcrosscov[0]: %g %g %g\r\n",tcrosscov[0][0],tcrosscov[0][1],tcrosscov[0][2]);
	fprintf(fptext,"tcrosscov[1]: %g %g %g\r\n",tcrosscov[1][0],tcrosscov[1][1],tcrosscov[1][2]);
	fprintf(fptext,"tcrosscov[2]: %g %g %g\r\n",tcrosscov[2][0],tcrosscov[2][1],tcrosscov[2][2]);
	return;
}
