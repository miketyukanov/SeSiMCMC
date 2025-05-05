/****************************************************************************\
SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001
$Id: StupidTest.cpp 1014 2009-03-01 16:50:36Z favorov $
\****************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

#include "../Sequences.hpp"
#include "../MarkovChainState.hpp"
#include "../KullbakCounter.hpp"

extern "C"
{
    #include "../Random.h"
}

unsigned int fake_data_positions[20]={84,81,78,75,72,69,66,63,60,50,
													47,44,41,38,35,32,29,26,23,20};

void createFakeExperimentalBunch(SequencesPile &sp,unsigned int be_quiet=0,
																	ostream * log_file_ptr=&cerr,
																	double noise_prob=0);

void createFakeExperimentalGappedBunch(SequencesPile &sp,unsigned int be_quiet=0,
																	ostream * log_file_ptr=&cerr,
																	double noise_prob=0);


int main(int argc, char ** argv)
{
/*****************************************************************************/
/*here, we put the defaults for all things, except those which depend on data*/ 
/*****************************************************************************/
	SequencesPile sp;
	rinit(1248312,2732596);
	createFakeExperimentalBunch(sp);
	MarkovChainState the_state (sp.size());
	vector<double> background(4);
	background[0]=
		background[1]=
		background[2]=
		background[3]=100;
	SymbolsCounter * Symbols = new AtgcSymbolsCounter(sp,the_state,1,background);
	KullbakCounter * Kullbak=new KullbakCounter(*Symbols);
	double long uuu;
	while(1)
	{
		uuu=Kullbak->EntropyThreshold();
	};
	delete Kullbak;
	
}

void createFakeExperimentalBunch(SequencesPile &sp,
																	unsigned int be_quiet,
																		ostream * log_file_ptr,
																		double noise_prob
																		)
//*sp is already created
{
	//20 sequnces, 105 random letters of ATGC (1234),
	// pattern ATGGCCATAA with random change in one letter with probability prob
	// is included in positions: 95,90,85......0 of corresponding pattern.
	unsigned int seq_no=20;
	unsigned int seq_length=120;
	unsigned int p_leng=10;
	unsigned short tpattern[10]={1,2,3,3,4,4,1,2,3,1};
	unsigned short cpattern[10]={2,4,1,2,3,3,4,4,1,2};
	unsigned short *pattern;
	if (!be_quiet) *log_file_ptr<<"Positions of motifs in test bunch are:\n";
	for (unsigned int i=0;i<seq_no;i++)
	{
		unsigned short * a=(unsigned short*)calloc(seq_length+1,sizeof(unsigned short));
		a[seq_length]=0;
		pattern=(i==8)?cpattern:tpattern; //9 has complemeted seq
//		pattern=(i==8)?tpattern:tpattern; //9 has not complemeted seq
		for (unsigned int j=0;j<seq_length;j++) a[j]=1+(int)floor(4*uni());
		//random seq is ready. Adding pattern
		if (i!=4) //omitting 5-th sequence
		{
			if (uni()<noise_prob) //with noise or not
			{
				unsigned int ch_pos=(int)floor(p_leng*uni());
				unsigned short ch_let=1+(int)floor(4*uni());
				unsigned short old_let=pattern[ch_pos];
				pattern[ch_pos]=ch_let;
				//we have change random letter of the pattern,
				//saved the original to old_let
				for (unsigned int j=0;j<p_leng;j++) a[fake_data_positions[i]+j]=
					pattern[j];
				pattern[ch_pos]=old_let; //restore
			}
			else
			{
				for (unsigned int j=0;j<p_leng;j++) a[fake_data_positions[i]+j]=
					pattern[j];
			}
			if (!be_quiet) *log_file_ptr<<fake_data_positions[i]<<" ";
		}
		else
		{
			for (unsigned int j=0;j<seq_length;j++) a[j]=j%4+1;
			if (!be_quiet) *log_file_ptr<<fake_data_positions[i]<<"** ";
		}
		string name="Fake|seq";
		char str_num[5];
		sprintf(str_num,"%02i",i);
		name.append(str_num);
		//the sequence is prepared.
		sp.add(a,seq_length,name);
	}
	if (!be_quiet) *log_file_ptr<<endl;
}
