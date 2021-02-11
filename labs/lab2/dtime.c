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

//=======================timer stuff ===============
#include <stdio.h>
#if 0 // def WIN32
#include <windows.h>
#include <tchar.h>
#endif
#if 1 // def _LINUX 
#include <sys/time.h>
#endif

// Timing Prototypes: Include where needed
void  dtimeInit();
double dtimeGet();
double dtimeSince(double t1,char *str);
double dtimeElapse(double t1);

// Timing Definitions: Included once.
static double Global_Tickrate = 0.0;
#define DTIMEINITMSG _T("You need to call dtimeInit(); once on Windows")
#define DTIMEINITERR _T("Coding Error")

double dtimeGet()
{
	double t;
#if 0 //def WIN32
	LARGE_INTEGER tvalue;
	QueryPerformanceCounter(&tvalue);   
	t = (double)(tvalue.QuadPart);
#endif
#if 1 //def _LINUX
	struct timeval tv1;
	gettimeofday(&tv1,NULL); 
	t = (tv1.tv_sec) + 1.0e-6*tv1.tv_usec;
#endif
	return(t);
}

double dtimeSince(double t1,char *str) 
{
	double t2;
	double telapsed;
#if 0 //def WIN32
	LARGE_INTEGER tvalue;
	QueryPerformanceCounter(&tvalue);   
	t2 = (double)(tvalue.QuadPart);
	if( Global_Tickrate > 0.0 ) { telapsed = (t2-t1)/Global_Tickrate; }
	else  
	{ 
		telapsed = -1.0;
		MessageBox(NULL,DTIMEINITMSG,DTIMEINITERR,MB_OK);
	}
#endif
#if 1 //def _LINUX
	struct timeval tv2;
	gettimeofday(&tv2,NULL); 
	t2 = (tv2.tv_sec) + 1.0e-6*tv2.tv_usec;
	telapsed = t2 - t1;
#endif
	printf("%.5g secs <-- Elapsed Time for: '%s'\r\n",telapsed,str);
	fflush(stdout);
	return( telapsed );
}

double dtimeElapse(double t1) 
{
	double t2;
	double telapsed;
#if 0 //def WIN32
	LARGE_INTEGER tvalue;
	QueryPerformanceCounter(&tvalue);   
	t2 = (double)(tvalue.QuadPart);
	if( Global_Tickrate > 0.0 ) { telapsed = (t2-t1)/Global_Tickrate; }
	else  
	{ 
		telapsed = -1.0;
		MessageBox(NULL,DTIMEINITMSG,DTIMEINITERR,MB_OK);
	}
#endif
#if 1 //def _LINUX
	struct timeval tv2;
	gettimeofday(&tv2,NULL); 
	t2 = (tv2.tv_sec) + 1.0e-6*tv2.tv_usec;
	telapsed = t2 - t1;
#endif
	return( telapsed );
}


void dtimeInit()
{
#if 0 //def WIN32
	LARGE_INTEGER cpufreq;
	QueryPerformanceFrequency(&cpufreq); 
	Global_Tickrate = (double)(cpufreq.QuadPart);
#endif
	return;
}
//===================================== end of timer stuff
