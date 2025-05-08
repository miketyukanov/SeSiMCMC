/****************************************************************************\
SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2021
$Id$
\****************************************************************************/

#include <stdlib.h>
#include <math.h>

using namespace std;

#include "Atgc.hpp"
#include "SymbolsCounter.hpp"
#include "Sequences.hpp"
//****************************************************************************
//we allocate all counters with letters_in_alphabet+1
//to use 1-based letters and to reserve 0 for control
//purposes.
//****************************************************************************
SymbolsCounter::PWM::PWM(const SymbolsCounter &sc):
	the_PWM((sc.letters_in_alphabet+1)*sc.pattern_length),
	letters_in_alphabet(sc.letters_in_alphabet),
	pattern_length(sc.pattern_length)
{
	for (unsigned int m=0;m<pattern_length;m++)
		for (unsigned short n=1;n<=letters_in_alphabet;n++)
			the_PWM[m*letters_in_alphabet+n]=foreground_probability(m,n);
}

SymbolsCounter::PWM & SymbolsCounter::PWM::operator=(const SymbolsCounter &sc)
{
	if (pattern_length!=sc.pattern_length)
	{
		the_PWM.resize((sc.letters_in_alphabet+1)*sc.pattern_length);
		pattern_length=sc.pattern_length;
	}
	for (unsigned int m=0;m<pattern_length;m++)
		for (unsigned short n=1;n<=letters_in_alphabet;n++)
			the_PWM[m*letters_in_alphabet+n]=sc.foreground_probability(m,n);
	return *this;
}

double SymbolsCounter::PWM::foreground_probability(unsigned int m,unsigned int n)
			const
{
	return the_PWM[m*letters_in_alphabet+n];
}

//****************************************************************************
//we allocate all counters with letters_in_alphabet+1
//to use 1-based letters and to reserve 0 for control
//purposes.
//****************************************************************************
SymbolsCounter::SymbolsCounter
(
	const SequencesPile & sp,
	unsigned int letters_in_alphabet_par,
	unsigned short if_common_background_par,
	const vector<double> & background
):
	Sequences(sp),
	letters_in_alphabet(letters_in_alphabet_par),
	pattern_length(0),
	spacer_5_end(0),
	spacer_3_end(0),
	is_spaced(0),
	if_common_background(if_common_background_par),
	total_letters(0),
	current_total_letters(0),
	motifs_count(0),
	letters(letters_in_alphabet+1),
	letters_per_sequence(letters_in_alphabet+1),
	b(letters_in_alphabet+1),
	c(), //c is resized when pattern_length changed
	c5prev(letters_in_alphabet+1),
	c5prev_second(letters_in_alphabet+1),
	s(letters_in_alphabet+1),
	s5prev(letters_in_alphabet+1)
{
	fill(letters.begin(),letters.end(),0);

	vector<unsigned int> counters(letters);  //sequence-wise counter

	if(!Sequences.size())
        throw DumbException("Trying to start counts from empty data set.\n");

    if (background.size()!=0 && background.size()!=4)
        throw DumbException("The background is not 0 and not 4.\n");

    if (background.size()==4) if_common_background=1;

	for (SequencesPile::const_iterator seq=sp.begin();seq<sp.end();seq++)
	{
		total_letters+=seq->size();
		fill(counters.begin(),counters.end(),0); //it is for current sequence only
		for (vector<unsigned short>::const_iterator symb=seq->begin();
				symb<seq->end();symb++)
		{
			letters[*symb]++;
			counters[*symb]++;
		}
		total_letters-=counters[0]; //we do not count mask
		if (!if_common_background)
			for (unsigned short letter=1;letter<=letters_in_alphabet;letter++)
				letters_per_sequence[letter].push_back(counters[letter]);
		//we never need it if background is common
	}

	if (background.size()==4)
	{
		double background_sum=0;
		for (unsigned short letter=1;letter<=letters_in_alphabet;letter++)
			background_sum+=background[letter-1];
		for (unsigned short letter=1;letter<=letters_in_alphabet;letter++)
			s[letter]=(int)rint((total_letters*background[letter-1])/background_sum);
		total_letters=0;
		for (unsigned short letter=1;letter<=letters_in_alphabet;letter++)
			total_letters+=s[letter];
	}

	else //background.size()!=4
	{
		for (unsigned short letter=1;letter<=letters_in_alphabet;letter++)
			s[letter]=letters[letter];
	}
	//pseudocounts
	//B=(B>4.)?B:4.;
	B=sqrt((double)Sequences.size());
	for(unsigned short i=1;i<=letters_in_alphabet;i++)
		b[i]=B*s[i]/total_letters;

	current_total_letters=total_letters;


};

double SymbolsCounter::restore_pseudocounts_sum ()
{
	return change_pseudocounts_sum((double)Sequences.size());
};


double SymbolsCounter::change_pseudocounts_sum (double sum)
{
	double old_sum=B;
	B=sum;
	for(unsigned short i=1;i<=letters_in_alphabet;i++)
		b[i]=B*s[i]/total_letters;
	return old_sum;
};

void SymbolsCounter::calculate (const MarkovChainState & theState)
	//gets the pattern_count and nonpattern_count values from the full dataset
	//to get the pattern_count and nonpattern_count for the dataset
	//excluding current sequence use exclude_sequence.
	//Backwards, use include_sequence.
{
	//we suppose that pattern length is critically less than nonpattern,
	//so it is easier to count the nonpattern counters only, taking pattern
	//as differences.
	if (!if_common_background)
		for (unsigned short letter=1;letter<=letters_in_alphabet;letter++)
			s[letter]=letters[letter];
	//now, s has all statistics exclude current_sequence

	motifs_count=0;
	//we init it

	pattern_length=theState.pattern_length();
	spacer_5_end=theState.spacer_5_end();
	spacer_3_end=theState.spacer_3_end();
	is_spaced=theState.is_spaced();
	informative_part_length=
			pattern_length-(is_spaced?(spacer_3_end-spacer_5_end+1):0);
	//pattern or gap length can be changed
	c.resize(pattern_length*(letters_in_alphabet+1),0);
	c.clear();
	//We believe that it is easier to reallocate it than to wipe

	if (!if_common_background)
		for (unsigned int p=0;p<Sequences.size();p++)//sequence - by sequence
		{
			if (!theState.motif_present[p]) continue; 	//it do not
																									//take part in pattern
																									//counters
			vector<unsigned short>::const_iterator ptr=
					Sequences[p].begin()+theState.positions[p];
			//ptr looks at first pattern letter
			for(unsigned int i=0;i<pattern_length;i++)
			{
				if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
				s[ptr[i]]--; //nonpattern_count decreases
				c[i*letters_in_alphabet+ptr[i]]++;//pattern counter increases
			}
			motifs_count++; //we have added the motif
		}
	else
		for (unsigned int p=0;p<Sequences.size();p++)//sequence - by sequence
		{
			if (!theState.motif_present[p]) continue; 	//it do not
																									//take part in pattern
																									//counters
			vector<unsigned short>::const_iterator ptr=
					Sequences[p].begin()+theState.positions[p];
			//ptr looks at first pattern letter
			for(unsigned int i=0;i<pattern_length;i++)
			{
				if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
				c[i*letters_in_alphabet+ptr[i]]++;//pattern counter increases
			}
			motifs_count++; //we have added the motif
		}

}

void SymbolsCounter::exclude_sequence
		(
			const MarkovChainState & theState,
			unsigned int sequence
		)
//All we do is removal of sequence from the sum.
//Works only if the pattern_length is the same for
//the state and the SymbolsCounter.
//We remove the sequence gain from the counter
//as it was done in accordance with theState.
//If it is untruth, the result can be stange
//
{
	if (!if_common_background)
		for (unsigned short letter=1;letter<=letters_in_alphabet;letter++)
			s[letter]=s[letter]-letters_per_sequence[letter][sequence];
	//now, s has all nonpattern statistics of prev. step
	//minus all sequence

	//remove current_sequence from pattern stats, too, rebuilding
	//the nonpattern stats
	//
	current_total_letters=current_total_letters-Sequences[sequence].size();
	//we have removed the sequence from current total letters count
	if (theState.motif_present[sequence])
		//it do take part in pattern counters
	{
		vector<unsigned short>::const_iterator ptr=
				Sequences[sequence].begin()+theState.positions[sequence];
		//ptr looks at first pattern letter
		if (!if_common_background)
			for(unsigned int i=0;i<pattern_length;i++)
			{
				if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
				s[ptr[i]]++;//nonpattern_count re-increases, we were not to
											//remove the letter from nonpattern stats
				c[i*letters_in_alphabet+ptr[i]]--;
											//pattern counter decreases
			}
		else
			for(unsigned int i=0;i<pattern_length;i++)
			{
				if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
				c[i*letters_in_alphabet+ptr[i]]--;
											//pattern counter decreases
			}
		motifs_count--; //we have removed the motif
	}
}

void SymbolsCounter::include_sequence
		(
			const MarkovChainState & theState,
			unsigned int sequence
		)
//All we do is addition of current_sequence from the sum in
//accordance with theState.
//Works only if the pattern_length is the same for
//the state and the SymbolsCounter.
//
{
	if (!if_common_background)
		for (unsigned short letter=1;letter<=letters_in_alphabet;letter++)
			s[letter]=s[letter]+letters_per_sequence[letter][sequence];
	//now, s has all nonpattern statistics of prev. step
	//minus all current_sequenceplus all previous_current_sequence
	//
	//remove current_sequence from pattern stats, too, rebuilding
	//the nonpattern stats
	//previous currernt sequence pattern adding to the counters
	current_total_letters=current_total_letters+Sequences[sequence].size();
	//we have added back the sequence from current total letters count
	if (theState.motif_present[sequence])
		//it do take part in pattern counters
	{
		vector<unsigned short>::const_iterator ptr=
				Sequences[sequence].begin()+theState.positions[sequence];
		//ptr looks at first pattern letter
		if (!if_common_background)
			for(unsigned int i=0;i<pattern_length;i++)
			{
				if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
				s[ptr[i]]--;//nonpattern_count decreases, we were not to
											//add the letter from nonpattern stats
				c[i*letters_in_alphabet+ptr[i]]++;
											//pattern counter increases
			}
		else
			for(unsigned int i=0;i<pattern_length;i++)
			{
				if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
				c[i*letters_in_alphabet+ptr[i]]++;
											//pattern counter increases
			}
		motifs_count++; //we have added the motif
	}
}


void SymbolsCounter::calculate_incrementally_after_3_end_shift
							(const MarkovChainState & theState)
//gets the pattern_count and nonpattern_count values from the full dataset
//after a shift of all motifs to right for 1 position.
//theState is the state after the shift.
//Of course, prevoius action of the SymbolsCounter object was
//calculate_complete or calculate_incrementally_after_3_end_shift
//
//Then, the state is shifted as a solid bunch, and now we are to
//reflect it in the counters.
{
	//we suppose that pattern length is critically less than nonpattern,
	//so it is easier to count the nonpattern counters only, taking pattern
	//as differences.
	//now, s has all statistics
	//s s already counted
	// ......****...****....
	// .......****...****....
	//cout<<"###### ("<<is_spaced<<" "<<spacer_5_end<<" "<<spacer_3_end<<")"<<endl<<flush;
	for(unsigned int j=1;j<=letters_in_alphabet;j++)
	{
		c5prev[j]=c[j];
		if (is_spaced)
			c5prev_second[j]=c[(spacer_3_end+1)*letters_in_alphabet+j];
		//it is, the next after spacer_3_end, for it will fall in gap.

	}
	if (!if_common_background)
		copy(s.begin(),s.end(),s5prev.begin());
	//we will use the *prev when calculate entropy
	for(unsigned int j=1;j<=letters_in_alphabet;j++)
		for(unsigned int i=0;i<pattern_length;i++)
			c[i*letters_in_alphabet+j]=c[(i+1)*letters_in_alphabet+j];
	//pattern counters shifts left (to 5')
	for(unsigned int j=1;j<=letters_in_alphabet;j++)
	{
		c[(pattern_length)*letters_in_alphabet+j]=0;
	//fills the head (3') pattern position with zeroes,
	//it will be filled later
		if (is_spaced)
		{
			c[(spacer_3_end)*letters_in_alphabet+j]=0; //now, it is gap position
			c[(spacer_5_end-1)*letters_in_alphabet+j]=0;
			//it is head of left part; it comes from spaced area,
			//and we do not know what is written here.
		}
	}



	if (!if_common_background)
		for (unsigned int p=0;p<Sequences.size();p++)//sequence - by sequence
		{
			if (!theState.motif_present[p]) continue; 	//it do not
																									//take part in pattern
																									//counters
			vector<unsigned short>::const_iterator ptr=
					Sequences[p].begin()+theState.positions[p]-1;
			//ptr looks at letter before first pattern letter an
			//we remove it from the statistics - the patterns shifted away
			//the pattern positions goes from 5' to 3' and
			//new 3' end will be in pattern_length from ptr
			//and ptr itself shows to the letter which has
			//disappered from pattern space
			unsigned short dead_letter=ptr[0];
			unsigned short new_letter=ptr[pattern_length];
			s[dead_letter]++; //nonpattern_count re-increases
			//now, we count the rightmost pattern position
			s[new_letter]--; //nonpattern_count decreases
			c[(pattern_length-1)*letters_in_alphabet+new_letter]++;
			//pattern counter increases
			if (is_spaced)
			{
				unsigned short dead_letter=
						*(Sequences[p].begin()+theState.positions[p]+spacer_3_end);
				//the position that now is 3-end of gap previously was
				//tail of right part.
				unsigned short new_letter=
						*(Sequences[p].begin()+theState.positions[p]+spacer_5_end-1);
				s[dead_letter]++; //nonpattern_count re-increases
				//now, we count the rightmost pattern position
				s[new_letter]--; //nonpattern_count decreases
				c[(spacer_5_end-1)*letters_in_alphabet+new_letter]++;
				//pattern counter increases
			}
		}
	else //if_common_background
		for (unsigned int p=0;p<Sequences.size();p++)//sequence - by sequence
		{
			if (!theState.motif_present[p]) continue; 	//it do not
																									//take part in pattern
																									//counters
			vector<unsigned short>::const_iterator ptr=
					Sequences[p].begin()+theState.positions[p]-1;
			//ptr looks at letter before first pattern letter an
			//we remove it from the statistics - the patterns shifted away
			//the pattern positions goes from 5' to 3' and
			//new 3' end will be in pattern_length from ptr
			//and ptr itself shows to the letter which has
			//disappered from pattern space
			unsigned short new_letter=ptr[pattern_length];
			c[(pattern_length-1)*letters_in_alphabet+new_letter]++;
			//pattern counter increases
			if (is_spaced)
			{
				unsigned short new_letter=
						*(Sequences[p].begin()+theState.positions[p]+spacer_5_end-1);
				c[(spacer_5_end-1)*letters_in_alphabet+new_letter]++;
				//pattern counter increases
			}
		}

}

double SymbolsCounter::foreground_probability(unsigned int m,unsigned int n)
			const
//n-th letter in m-th motif position probability
//position is 0-based; letter is 1-based.
{
	if (is_spaced && m>=spacer_5_end && m<=spacer_3_end)
		return background_probability(n);
	return (((double)pattern_count(m,n)+pseudocount(n))/
	(motifs()+pseudocounts_sum()));
}

double SymbolsCounter::foreground_weight(unsigned int m,unsigned int n)
			const
//n-th letter in m-th motif position weight
//position is 0-based; letter is 1-based.
{
	if (is_spaced && m>=spacer_5_end && m<=spacer_3_end)
		return background_probability(n);
	return ((double)pattern_count(m,n)/motifs());
}

double SymbolsCounter::background_probability(unsigned int n) const
//n-th letter background probability
//position is 0-based; letter is 1-based.
{
	if (!if_common_background)
		return ((nonpattern_count(n)+pseudocount(n))/
			(total()+pseudocounts_sum()-
			informative_part_length*motifs()));
	else //common_background
		return ((/*letters[n]*/nonpattern_count(n)+pseudocount(n))/
			(total()+pseudocounts_sum()));
}

double long
SymbolsCounter::NegativeEntropy()
//F in Lawrence93
{
//	cout<<"\n__________________________________start...\n";
//	cout<<*this<<endl;
	double long PEntropy=0;
	for (unsigned int i=0;i<pattern_length;i++)
	{
//		double long P1=PEntropy;
		if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
		for (unsigned short j=1;j<=letters_in_alphabet;j++)
			PEntropy+=pattern_count(i,j)*
							log_2(foreground_probability(i,j)/background_probability(j));
//		cout<<(PEntropy-P1)<<"  ";
	}
//	cout<<"\n__________________________________fin...("<<PEntropy<<")\n";
	return PEntropy;
}



double long
SymbolsCounter::NegativeEntropy_change_after_3_end_shift()
{
	//the change contains of 3 parts -
	//old 5'-enddisappears, new 3'-end appears
	//and the middle part changes because of
	//background counts change.
	double long PEntropyChange=0.;
	for (unsigned short j=1;j<=letters_in_alphabet;j++)
	{
		double q,p,p_old;

		q=
				(c5prev[j]+pseudocount(j))/
				(motifs()+pseudocounts_sum());
		//actually, p_old==p if if_common_background
		p=
				(nonpattern_count(j)+pseudocount(j))/
				(total()+pseudocounts_sum()-
				 pattern_length*motifs());

		if (!if_common_background)
			p_old=
					(s5prev[j]+pseudocount(j))/
					(total()+pseudocounts_sum()-
					 pattern_length*motifs());
		else p_old=p;

		PEntropyChange-=c5prev[j]*log_2(q/p_old);

		if (is_spaced)
		{
			q=
					(c5prev_second[j]+pseudocount(j))/
					(motifs()+pseudocounts_sum());
			PEntropyChange-=c5prev_second[j]*log_2(q/p_old);
		}

		q=
				(pattern_count(pattern_length-1,j)+pseudocount(j))/
				(motifs()+pseudocounts_sum());

		PEntropyChange+=pattern_count(pattern_length-1,j)*log_2(q/p);

		if (is_spaced)
		{
			q=
					(pattern_count(spacer_5_end-1,j)+pseudocount(j))/
					(motifs()+pseudocounts_sum());

			PEntropyChange+=pattern_count(spacer_5_end-1,j)*log_2(q/p);
		}

		if (!if_common_background)
			for (unsigned int i=0;i<pattern_length-1;i++)
			{
				if (is_spaced && i==spacer_5_end-1) i=spacer_3_end+1;
				PEntropyChange-=pattern_count(i,j)*log_2(p/p_old);
			}
	}
	return PEntropyChange;
}

/****************************************************************************
 ****************************************************************************/
SymmetricSymbolsCounter::SymmetricSymbolsCounter
	(
		const SequencesPile & sp,
		unsigned int letters_in_alphabet_par,
		unsigned short if_common_background_par=0,
		const vector<double> & background=*new vector<double>
	):
	SymbolsCounter(sp,letters_in_alphabet_par,if_common_background_par)
{
	NegativeEntropyValue=0;
};

double SymmetricSymbolsCounter::foreground_probability(unsigned int m,unsigned int n)
			const
//n-th letter in m-th motif position probability
//position is 0-based; letter is 1-based.
{
	if (is_spaced && m>=spacer_5_end && m<=spacer_3_end)
		return background_probability(n);
	unsigned short nc=Atgc::complement((short unsigned)n);
	return ((double)pattern_count(m,n)+(double)pattern_count(pattern_length-m-1,nc)+
	 pseudocount(n)+pseudocount(nc))/
	(2*(motifs()+pseudocounts_sum()));
}

double SymmetricSymbolsCounter::foreground_weight(unsigned int m,unsigned int n)
			const
//n-th letter in m-th motif position weight
//position is 0-based; letter is 1-based.
{
	if (is_spaced && m>=spacer_5_end && m<=spacer_3_end)
		return background_probability(n);
	unsigned short nc=Atgc::complement((short unsigned)n);
	return ((double)pattern_count(m,n)+(double)pattern_count(pattern_length-m-1,nc))/
	(2*(motifs()));
}

double long
SymmetricSymbolsCounter::NegativeEntropy()
//F in Lawrence93
//it is for palindroms, it do not have
//incremental version.
//the idea is that still we have usual
//SymbolsCounter counter matrix,
//the actual foreground probability is
//the mean for q(letter,pos) and q(complementary_letter,symmetric_pos)

{
	double long PEntropy=0;
	for (unsigned int i=0;i<pattern_length;i++)
	{
		if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
		for (unsigned short j=1;j<=letters_in_alphabet;j++)
			PEntropy+=pattern_count(i,j)*
					log_2(foreground_probability(i,j)/background_probability(j));
	}
	NegativeEntropyValue=PEntropy;
	return PEntropy;
}

double long
SymmetricSymbolsCounter::NegativeEntropy_change_after_3_end_shift()
{
	long double pe=NegativeEntropyValue;
	return NegativeEntropy()-pe;
}
/****************************************************************************
 ****************************************************************************/
DoubletSymbolsCounter::DoubletSymbolsCounter
	(
		const SequencesPile & sp,
		unsigned int letters_in_alphabet_par,
		unsigned short if_common_background_par=0,
		const vector<double> & background=*new vector<double>
	):
	SymbolsCounter(sp,letters_in_alphabet_par,if_common_background_par,background)
{
	NegativeEntropyValue=0;
};

double DoubletSymbolsCounter::foreground_probability(unsigned int m,unsigned int n)
			const
//n-th letter in m-th motif position probability
//position is 0-based; letter is 1-based.
{
	unsigned int m_second=m;
	if (is_spaced)
	{
		if(m>=spacer_5_end && m<=spacer_3_end)
			return background_probability(n);
		m_second=(pattern_length-m>m)?spacer_3_end+1+m:m-spacer_3_end-1;
	}
	else
		m_second=(pattern_length-m>m)?m+((pattern_length+1)/2):m-((pattern_length+1)/2);

	return ((double)pattern_count(m,n)+(double)pattern_count(m_second,n)+
	 pseudocount(n)+pseudocount(n))/
	(2*(motifs()+pseudocounts_sum()));
}

double DoubletSymbolsCounter::foreground_weight(unsigned int m,unsigned int n)
			const
//n-th letter in m-th motif position weight
//position is 0-based; letter is 1-based.
{
	unsigned int m_second=m;
	if (is_spaced)
	{
		if(m>=spacer_5_end && m<=spacer_3_end)
			return background_probability(n);
		m_second=(pattern_length-m>m)?spacer_3_end+1+m:m-spacer_3_end-1;
	}
	else
		m_second=(pattern_length-m>m)?m+((pattern_length+1)/2):m-((pattern_length+1)/2);

	return ((double)pattern_count(m,n)+(double)pattern_count(m_second,n))/
	(2*(motifs()));
}

double long
DoubletSymbolsCounter::NegativeEntropy()
//F in Lawrence93
//it is for palindroms, it do not have
//incremental version.
//the idea is that still we have usual
//SymbolsCounter counter matrix,
//the actual foreground probability is
//the mean for q(letter,pos) and q(complementary_letter,symmetric_pos)

{
	double long PEntropy=0;
	for (unsigned int i=0;i<pattern_length;i++)
	{
		if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
		for (unsigned short j=1;j<=letters_in_alphabet;j++)
			PEntropy+=pattern_count(i,j)*
					log_2(foreground_probability(i,j)/background_probability(j));
	}
	NegativeEntropyValue=PEntropy;
	return PEntropy;
}

double long
DoubletSymbolsCounter::NegativeEntropy_change_after_3_end_shift()
{
	long double pe=NegativeEntropyValue;
	return NegativeEntropy()-pe;
}

/****************************************************************************
 ****************************************************************************/
DoubletAtgcSymbolsCounter::DoubletAtgcSymbolsCounter
	(
		const SequencesPile & sp,
		const MarkovChainState & mcs,
		unsigned short if_common_background_par=0,
		const vector<double> & background=*new vector<double>
	):
	AtgcSymbolsCounter(sp,mcs,if_common_background_par,background)
{
	NegativeEntropyValue=0;
};

double DoubletAtgcSymbolsCounter::foreground_probability(unsigned int m,unsigned int n)
			const
//n-th letter in m-th motif position probability
//position is 0-based; letter is 1-based.
{
	unsigned int m_second=m;
	if (is_spaced)
	{
		if(m>=spacer_5_end && m<=spacer_3_end)
			return background_probability(n);
		m_second=(pattern_length-m>m)?spacer_3_end+1+m:m-spacer_3_end-1;
	}
	else
		m_second=(pattern_length-m>m)?m+((pattern_length+1)/2):m-((pattern_length+1)/2);

	return ((double)pattern_count(m,n)+(double)pattern_count(m_second,n)+
	 pseudocount(n)+pseudocount(n))/
	(2*(motifs()+pseudocounts_sum()));
}

double DoubletAtgcSymbolsCounter::foreground_weight(unsigned int m,unsigned int n)
			const
//n-th letter in m-th motif position weight
//position is 0-based; letter is 1-based.
{
	unsigned int m_second=m;
	if (is_spaced)
	{
		if(m>=spacer_5_end && m<=spacer_3_end)
			return background_probability(n);
		m_second=(pattern_length-m>m)?spacer_3_end+1+m:m-spacer_3_end-1;
	}
	else
		m_second=(pattern_length-m>m)?m+((pattern_length+1)/2):m-((pattern_length+1)/2);

	return ((double)pattern_count(m,n)+(double)pattern_count(m_second,n))/
	(2*(motifs()));
}

double long
DoubletAtgcSymbolsCounter::NegativeEntropy()
//F in Lawrence93
//it is for palindroms, it do not have
//incremental version.
//the idea is that still we have usual
//SymbolsCounter counter matrix,
//the actual foreground probability is
//the mean for q(letter,pos) and q(complementary_letter,symmetric_pos)

{
	double long PEntropy=0;
	for (unsigned int i=0;i<pattern_length;i++)
	{
		if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
		for (unsigned short j=1;j<=letters_in_alphabet;j++)
			PEntropy+=pattern_count(i,j)*
					log_2(foreground_probability(i,j)/background_probability(j));
	}
	NegativeEntropyValue=PEntropy;
	return PEntropy;
}

double long
DoubletAtgcSymbolsCounter::NegativeEntropy_change_after_3_end_shift()
{
	long double pe=NegativeEntropyValue;
	return NegativeEntropy()-pe;
}


/****************************************************************************
 ****************************************************************************/

AtgcSymbolsCounter::AtgcSymbolsCounter
(
	const SequencesPile & sp,
	const MarkovChainState & mcs,
	unsigned short if_common_background_par=0,
	const vector<double> & background=*new vector<double>
):
SymbolsCounter(sp,4,if_common_background_par,background),
counted_as_complement(sp.size())
{
	fill(counted_as_complement.begin(),counted_as_complement.end(),0);
	for (unsigned int p=0;p<sp.size();p++)
		if (mcs.is_complement[p]) current_sequence_was_complemented(p);
}

void AtgcSymbolsCounter::current_sequence_was_complemented(unsigned int seq_no)
{
	counted_as_complement[seq_no]=!counted_as_complement[seq_no];
	if (if_common_background) return;
	unsigned int letter,complement_letter;
	letter=1;
	complement_letter=2;
	{
		letters[letter]=letters[letter]
				-letters_per_sequence[letter][seq_no]
				+letters_per_sequence[complement_letter][seq_no];
		letters[complement_letter]=letters[complement_letter]
				-letters_per_sequence[complement_letter][seq_no]
				+letters_per_sequence[letter][seq_no];
		swap(
			letters_per_sequence[letter][seq_no],
			letters_per_sequence[complement_letter][seq_no]
			);
	};
	letter=3;
	complement_letter=4;
	{
		letters[letter]=letters[letter]
				-letters_per_sequence[letter][seq_no]
				+letters_per_sequence[complement_letter][seq_no];
		letters[complement_letter]=letters[complement_letter]
				-letters_per_sequence[complement_letter][seq_no]
				+letters_per_sequence[letter][seq_no];
		swap(
			letters_per_sequence[letter][seq_no],
			letters_per_sequence[complement_letter][seq_no]
			);
	};
//The two {} are identical, but it makes no sense to define a
//class-level function to "swap data in r and pc for this sequence for
//these two letters and call it twice :)
	for(unsigned short i=1;i<=letters_in_alphabet;i++)
		b[i]=B*letters[i]/total_letters;
//pseudocounts
}

void AtgcSymbolsCounter::calculate (const MarkovChainState & theState)
{

	//we suppose that pattern length is critically less than nonpattern,
	//so it is easier to count the nonpattern counters only, taking pattern
	//as differences.
	for (unsigned int p=0;p<Sequences.size();p++)
		if(theState.is_complement[p]!=counted_as_complement[p])
			current_sequence_was_complemented(p);

	if (!if_common_background)
	{
		for (unsigned short letter=1;letter<=letters_in_alphabet;letter++)
			s[letter]=letters[letter];
	}
	//now, s has all statistics exclude current_sequence

	motifs_count=0; //we init it

	pattern_length=theState.pattern_length();
	spacer_5_end=theState.spacer_5_end();
	spacer_3_end=theState.spacer_3_end();
	is_spaced=theState.is_spaced();
	informative_part_length=
			pattern_length-(is_spaced?(spacer_3_end-spacer_5_end+1):0);

	//pattern length can be changed
	c.resize(pattern_length*(letters_in_alphabet+1),0);
	c.clear();
	//We believe that it is easier to reallocate it than to wipe

	if (!if_common_background)
		for (unsigned int p=0;p<Sequences.size();p++)//sequence - by sequence
		{
			if (!theState.motif_present[p]) continue; 	//it do not
																									//take part in pattern
																									//counters

			vector<unsigned short>::const_iterator ptr=
					Sequences[p].begin()+theState.positions[p];
			//ptr looks at first pattern letter
			if (theState.is_complement[p])
				for(unsigned int i=0;i<pattern_length;i++)
				{
					if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
					unsigned short letter=Atgc::complement(ptr[pattern_length-1-i]);
					s[letter]--; //nonpattern_count decreases
					c[i*letters_in_alphabet+letter]++;//pattern counter increases
				}
			else
				for(unsigned int i=0;i<pattern_length;i++)
				{
					if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
					unsigned short letter=ptr[i];

					s[letter]--; //nonpattern_count decreases

					c[i*letters_in_alphabet+letter]++;//pattern counter increases

				}


			motifs_count++; //we have added the motif
		}
	else
		for (unsigned int p=0;p<Sequences.size();p++)//sequence - by sequence
		{
			if (!theState.motif_present[p]) continue; 	//it do not
																									//take part in pattern
																									//counters
			vector<unsigned short>::const_iterator ptr=
					Sequences[p].begin()+theState.positions[p];
			//ptr looks at first pattern letter
			if (theState.is_complement[p])
				for(unsigned int i=0;i<pattern_length;i++)
				{
					if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
					unsigned short letter=Atgc::complement(ptr[pattern_length-1-i]);
					c[i*letters_in_alphabet+letter]++;//pattern counter increases
				}
			else
				for(unsigned int i=0;i<pattern_length;i++)
				{
					if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
					unsigned short letter=ptr[i];
					c[i*letters_in_alphabet+letter]++;//pattern counter increases
				}
			motifs_count++; //we have added the motif
		}

}

void AtgcSymbolsCounter::exclude_sequence
		(
			const MarkovChainState & theState,
			unsigned int sequence
		)
{
	if (!if_common_background)
		for (unsigned short letter=1;letter<=letters_in_alphabet;letter++)
			s[letter]=s[letter]-letters_per_sequence[letter][sequence];
	//now, s has all nonpattern statistics of prev. step
	//minus all current_sequenceplus all previous_current_sequence
	//
	//remove current_sequence from pattern stats, too, rebuilding
	//the nonpattern stats
	current_total_letters=current_total_letters-Sequences[sequence].size();
	//we have removed the sequence from current total letters count
	if (theState.motif_present[sequence])
		//it do take part in pattern counters
	{
		vector<unsigned short>::const_iterator ptr=
				Sequences[sequence].begin()+theState.positions[sequence];
		//ptr looks at first pattern letter
		if (!if_common_background)
		{
			if (theState.is_complement[sequence])
				for(unsigned int i=0;i<pattern_length;i++)
				{
					if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
					unsigned short letter=Atgc::complement(ptr[pattern_length-1-i]);
					s[letter]++; //nonpattern_count re-increases
					c[i*letters_in_alphabet+letter]--;//pattern counter re-decreases
				}
			else
				for(unsigned int i=0;i<pattern_length;i++)
				{
					if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
					unsigned short letter=ptr[i];
					s[letter]++; //nonpattern_count re-increases
					c[i*letters_in_alphabet+letter]--;//pattern counter re-decreases
				}
		}
		else //if_common_background
		{
			if (theState.is_complement[sequence])
				for(unsigned int i=0;i<pattern_length;i++)
				{
					if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
					unsigned short letter=Atgc::complement(ptr[pattern_length-1-i]);
					c[i*letters_in_alphabet+letter]--;//pattern counter re-decreases
				}
			else
				for(unsigned int i=0;i<pattern_length;i++)
				{
					if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
					unsigned short letter=ptr[i];
					c[i*letters_in_alphabet+letter]--;//pattern counter re-decreases
				}
		}

		motifs_count--; //we have removed the motif
	}

}
void AtgcSymbolsCounter::include_sequence
		(
			const MarkovChainState & theState,
			unsigned int sequence
		)
{
	if (!if_common_background)
		for (unsigned short letter=1;letter<=letters_in_alphabet;letter++)
			s[letter]=s[letter]+letters_per_sequence[letter][sequence];
	//now, s has all nonpattern statistics of prev. step
	//plus all sequence
	//
	//remove current_sequence from pattern stats, too, rebuilding
	//the nonpattern stats

	current_total_letters=current_total_letters+Sequences[sequence].size();
	//we have added back the sequence from current total letters count
	//previous currernt sequence pattern adding to the counters
	if (theState.motif_present[sequence])
		//it do take part in pattern counters
	{
		vector<unsigned short>::const_iterator ptr=
				Sequences[sequence].begin()+theState.positions[sequence];
		//ptr looks at first pattern letter
		if (!if_common_background)
		{
			if (theState.is_complement[sequence])
				for(unsigned int i=0;i<pattern_length;i++)
				{
					if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
					unsigned short letter=Atgc::complement(ptr[pattern_length-1-i]);
					s[letter]--; //nonpattern_count re-decreases
					c[i*letters_in_alphabet+letter]++;//pattern counter increases
				}
			else
				for(unsigned int i=0;i<pattern_length;i++)
				{
					if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
					unsigned short letter=ptr[i];
					s[letter]--; //nonpattern_count re-decreases
					c[i*letters_in_alphabet+letter]++;//pattern counter increases
				}
		}
		else //if_common_background
		{
			if (theState.is_complement[sequence])
				for(unsigned int i=0;i<pattern_length;i++)
				{
					if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
					unsigned short letter=Atgc::complement(ptr[pattern_length-1-i]);
					c[i*letters_in_alphabet+letter]++;//pattern counter increases
				}
			else
				for(unsigned int i=0;i<pattern_length;i++)
				{
					if (is_spaced && i==spacer_5_end) i=spacer_3_end+1;
					unsigned short letter=ptr[i];
					c[i*letters_in_alphabet+letter]++;//pattern counter increases
				}
		}
		motifs_count++; //we have added the motif
	}
}

void AtgcSymbolsCounter::calculate_incrementally_after_3_end_shift
							(const MarkovChainState & theState)
//gets the pattern_count and nonpattern_count values from the full dataset
//after a shift of all motifs to right for 1 position.
//theState is the state after the shift.
//Of course, prevoius action of the SymbolsCounter object was
//calculate_complete or calculate_complete_incrementally_after_right_shift
//
{
	//we suppose that pattern length is critically less than nonpattern,
	//so it is easier to count the nonpattern counters only, taking pattern
	//as differences.
	//now, s has all statistics
	//s s already counted

	for(unsigned int j=1;j<=letters_in_alphabet;j++)
	{
		c5prev[j]=c[j];
		if (is_spaced)
			c5prev_second[j]=c[(spacer_3_end+1)*letters_in_alphabet+j];
		//it is, the next after spacer_3_end, for it will fall in gap.
	}
	copy(s.begin(),s.end(),s5prev.begin());
	for(unsigned int j=1;j<=letters_in_alphabet;j++)
		for(unsigned int i=0;i<pattern_length;i++)
			c[i*letters_in_alphabet+j]=c[(i+1)*letters_in_alphabet+j];
	//pattern counters shifts left (to 5')
	for(unsigned int j=1;j<=letters_in_alphabet;j++)
	{
		c[(pattern_length-1)*letters_in_alphabet+j]=0;
	//fills the last (3') pattern position with zeroes,
	//it will be filled later
		if (is_spaced)
		{
			c[(spacer_3_end)*letters_in_alphabet+j]=0; //now, it is gap position
			c[(spacer_5_end-1)*letters_in_alphabet+j]=0;
			//it is head of left part; it comes from spaced area,
			//and we do not know what is written here.
		}
	}

	if (!if_common_background)
	{
		for (unsigned int p=0;p<Sequences.size();p++)//sequence - by sequence
		{
			if (!theState.motif_present[p]) continue; 	//it do not
																									//take part in pattern
																									//counters
			if (theState.is_complement[p])
			{
				//it is complement. So its 3'-shift was 5' in our space.
				vector<unsigned short>::const_iterator ptr=
						Sequences[p].begin()+theState.positions[p];
				//ptr looks at complemnt to last (3') pattern letter
				//it is the thing to be added to statistics
				//the previous 5'-end which was left by the pattern
				//is complement to prt[pattern_length]
				//we remove it from the statistics - the patterns shifted away
				unsigned short dead_letter=Atgc::complement(ptr[pattern_length]);
				unsigned short new_letter=Atgc::complement(ptr[0]);
				s[dead_letter]++; //nonpattern_count re-increases
				//now, we count the rightmost pattern position
				s[new_letter]--; //nonpattern_count decreases
				c[(pattern_length-1)*letters_in_alphabet+new_letter]++;
				//pattern counter increases
			}
			else
			{
				//it is noncomplement.
				vector<unsigned short>::const_iterator ptr=
						Sequences[p].begin()+theState.positions[p]-1;
				//ptr looks at letter before first pattern letter an
				//we remove it from the statistics - the patterns shifted away
				//the pattern positions goes from 5' to 3' and
				//new 3' end will be in pattern_length from ptr
				//and ptr itself shows to the letter which has
				//disappered from pattern space
				unsigned short dead_letter=ptr[0];
				unsigned short new_letter=ptr[pattern_length];
				s[dead_letter]++; //nonpattern_count re-increases
				//now, we count the rightmost pattern position
				s[new_letter]--; //nonpattern_count decreases
				c[(pattern_length-1)*letters_in_alphabet+new_letter]++;
				//pattern counter increases
			}
			if (is_spaced)
			{
				if (theState.is_complement[p])
				{
					//it is complement. So its 3'-shift was 5' in our space.
					unsigned short dead_letter=
							Atgc::complement
								(*(Sequences[p].begin()+theState.positions[p]+spacer_5_end));
					unsigned short new_letter=
							Atgc::complement
								(*(Sequences[p].begin()+theState.positions[p]+spacer_3_end+1));
					s[dead_letter]++; //nonpattern_count re-increase
					//now, we count the rightmost pattern position
					s[new_letter]--; //nonpattern_count decreases
					c[(spacer_5_end-1)*letters_in_alphabet+new_letter]++;
					//pattern counter increases
				}
				else
				{
					//it is noncomplement.
					unsigned short dead_letter=
							*(Sequences[p].begin()+theState.positions[p]+spacer_3_end);
					unsigned short new_letter=
							*(Sequences[p].begin()+theState.positions[p]+spacer_5_end-1);
					s[dead_letter]++; //nonpattern_count re-increase
					//now, we count the rightmost pattern position
					s[new_letter]--; //nonpattern_count decreases
					c[(spacer_5_end-1)*letters_in_alphabet+new_letter]++;
					//pattern counter increases
				}
			}
		}
	}
	else
	{
		for (unsigned int p=0;p<Sequences.size();p++)//sequence - by sequence
		{
			if (!theState.motif_present[p]) continue; 	//it do not
																									//take part in pattern
																									//counters
			if (theState.is_complement[p])
			{
				//it is complement. So its 3'-shift was 5' in our space.
				vector<unsigned short>::const_iterator ptr=
						Sequences[p].begin()+theState.positions[p];
				//ptr looks at complemnt to last (3') pattern letter
				//it is the thing to be added to statistics
				//the previous 5'-end which was left by the pattern
				//is complement to prt[pattern_length]
				//we remove it from the statistics - the patterns shifted away
				unsigned short new_letter=Atgc::complement(ptr[0]);
				c[(pattern_length-1)*letters_in_alphabet+new_letter]++;
				//pattern counter increases
			}
			else
			{
				//it is noncomplement.
				vector<unsigned short>::const_iterator ptr=
						Sequences[p].begin()+theState.positions[p]-1;
				//ptr looks at letter before first pattern letter an
				//we remove it from the statistics - the patterns shifted away
				//the pattern positions goes from 5' to 3' and
				//new 3' end will be in pattern_length from ptr
				//and ptr itself shows to the letter which has
				//disappered from pattern space
				unsigned short new_letter=ptr[pattern_length];
				c[(pattern_length-1)*letters_in_alphabet+new_letter]++;
				//pattern counter increases
			}
			if (is_spaced)
			{
				if (theState.is_complement[p])
				{
					//it is complement. So its 3'-shift was 5' in our space.
					unsigned short new_letter=
							Atgc::complement
								(*(Sequences[p].begin()+theState.positions[p]+spacer_3_end+1));
					c[(spacer_5_end-1)*letters_in_alphabet+new_letter]++;
					//pattern counter increases
				}
				else
				{
					//it is noncomplement.
					unsigned short new_letter=
							*(Sequences[p].begin()+theState.positions[p]+spacer_5_end-1);
					c[(spacer_5_end-1)*letters_in_alphabet+new_letter]++;
					//pattern counter increases
				}
			}
		}
	}

}


