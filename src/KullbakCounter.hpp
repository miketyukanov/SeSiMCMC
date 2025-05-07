/****************************************************************************\
SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2021
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For general describtion of the classes declared in the header, see headers.txt
$Id$
\****************************************************************************/

#ifndef _KULLBAK_COUNTER_HPP
#define _KULLBAK_COUNTER_HPP

#include <math.h>

#include <vector>
#include <iostream>
#include <iomanip>

#include "Exception.hpp"
#include "SymbolsCounter.hpp"

class KullbakCounter
{
    KullbakCounter(const KullbakCounter &) = delete;
    KullbakCounter &operator=(const KullbakCounter &) = delete;
    //no copy constructor, no copy operation
	double entropy_threshold;
	//maximal entropy for a "stable" position :)
	//Pay attention: it is entropy, not IC!
	const SymbolsCounter & Symbols;
	vector<unsigned short> if_strong_position; //0 or 1
	double long strong_pos_IC;
	unsigned int strong_pos_counter;
	void weaken_position(unsigned int pos);
public:

	KullbakCounter(const SymbolsCounter & sc);

	void weaken_outer_strong_islands(unsigned int multiplier);

	static double long EntropyThreshold (const vector<double> & bg_probabilities);
	//bg_probabilities has size letters_in_alphabet+1,
	//[0] is not used
	static double long EntropyThreshold (const SymbolsCounter & sc)
	{
		vector<double> bg(sc.letters_in_alphabet+1);
		for (unsigned int i=1;i<=sc.letters_in_alphabet;i++)
			bg[i]=sc.background_probability(i);
		return EntropyThreshold(bg);
	};
	double long EntropyThreshold ()
	{
		return EntropyThreshold(Symbols);
	};

	double long EntropyDistanceFrom(const SymbolsCounter::PWM &) const;
	//count distance from parameter to instance distribution

	const vector<unsigned short> & if_a_strong_position; //0 or 1

	const double long & strong_positions_IC;
	const unsigned int & strong_positions_counter;
 	const unsigned int & informative_part_length;

};

inline
KullbakCounter::KullbakCounter(const SymbolsCounter & sc):
					entropy_threshold(0),Symbols(sc),
					if_strong_position(sc.pattern_length),
					strong_pos_IC(0),
					strong_pos_counter(0),
					if_a_strong_position(if_strong_position),
					strong_positions_IC(strong_pos_IC),
					strong_positions_counter(strong_pos_counter),
					informative_part_length(Symbols.informative_part_length)
{
	entropy_threshold=EntropyThreshold();
	for (unsigned int i=0;i<Symbols.pattern_length;i++) //what position are to be capped
	{
		double long entropy=0;
		if (Symbols.is_spaced && i==Symbols.spacer_5_end) i=Symbols.spacer_3_end+1;
		for (unsigned short j=1;j<=Symbols.letters_in_alphabet;j++)
			if (Symbols.foreground_weight(i,j)>EPSILON)
			{
					entropy-=Symbols.foreground_probability(i,j)*
					log_2(Symbols.foreground_probability(i,j));
			}
		strong_pos_counter+=if_strong_position[i]=(entropy<entropy_threshold);
		if (if_strong_position[i])
		{
			for (unsigned short j=1;j<=Symbols.letters_in_alphabet;j++)
				strong_pos_IC+=Symbols.foreground_weight(i,j)*
								log_2(Symbols.foreground_probability(i,j)/Symbols.background_probability(j));
		};
	};
};

inline
void
KullbakCounter::weaken_position(unsigned int pos)
{
	if_strong_position[pos]=0;
	strong_pos_counter--;
	for (unsigned short j=1;j<=Symbols.letters_in_alphabet;j++)
		strong_pos_IC-=Symbols.foreground_weight(pos,j)*
						log_2(Symbols.foreground_probability(pos,j)/Symbols.background_probability(j));

}

inline
void
KullbakCounter::weaken_outer_strong_islands(unsigned int multiplier)
{
	//if an outer island of strong positions is isolated by more than
	//island_length*multiplier of weak positions, it is wiped from strong positions
	unsigned int
		current_island_start=0,
		current_island_length=0,
		putative_island_start=0,
		putative_island_length=0,
		wing_finished=0;
	for (unsigned int pos=0; !wing_finished ;pos++) //left wing
	{
		if
		(
			(Symbols.is_spaced && pos==Symbols.spacer_5_end)
			||
			(pos+1>Symbols.pattern_length/2)
		)
			wing_finished=1;
		else //we analyse the current position
		{
			if (current_island_length==0) //we are looking for a new island
			{
				if (if_strong_position[pos]) //we start a new island
				{
					current_island_start=pos;
					current_island_length=1;
				}
			}
			else //we have a current island
			{
				if (if_strong_position[pos])
					current_island_length++; //we elongate it
				else //we finish it, now it is putative
				{
					putative_island_length=current_island_length;
					putative_island_start=current_island_start;
					current_island_length=0; //we do not have a current island
				}
			}
		}

		if (putative_island_length && (current_island_length || wing_finished))
			//we have a putative island and we have started a current, i.e. meet a strong position
			//or the wing has finished, so we are to decide what to do with the
			//memorised putative island
		{
			if (
						putative_island_length*multiplier<
						(pos-(putative_island_start+putative_island_length)-wing_finished)
					)
				//the stupid putative island is multiplier times less than the gap
				//between it and the next island; we remove it
				for (unsigned int weaken_pos=putative_island_start;
							weaken_pos<putative_island_start+putative_island_length;
							weaken_pos++) weaken_position(weaken_pos);
			putative_island_length=0; //we do not have a putative_island any more
		}
	}

	current_island_start=0,
	current_island_length=0,
	putative_island_start=0,
	putative_island_length=0,
	wing_finished=0;
	for (unsigned int pos=Symbols.pattern_length-1; !wing_finished ;pos--) //right wing
	{
		if
		(
			(Symbols.is_spaced && pos==Symbols.spacer_3_end)
			||
			(pos-1<Symbols.pattern_length/2)
		)
			wing_finished=1;
		else //we analyse the current position
		{
			if (current_island_length==0) //we are looking for a new island
			{
				if (if_strong_position[pos]) //we start a new island
				{
					current_island_start=pos;
					current_island_length=1;
				}
			}
			else //we have a current island
			{
				if (if_strong_position[pos])
					current_island_length++; //we elongate it
				else //we finish it, now it is putative
				{
					putative_island_length=current_island_length;
					putative_island_start=current_island_start;
					current_island_length=0; //we do not have a current island
				}
			}
		}

		if (putative_island_length && (current_island_length || wing_finished))
			//we have a putative island and we have started a current, i.e. meet a strong position
			//or the wing has finished, so we are to decide what to do with the
			//memorised putative island
		{
			if (
						putative_island_length*multiplier<
						(putative_island_start-putative_island_length-pos+wing_finished)
					)
				//the stupid putative island is multiplier times less than the gap
				//between it and the next island; we remove it
				for (unsigned int weaken_pos=putative_island_start;
							weaken_pos>putative_island_start-putative_island_length;
							weaken_pos--) weaken_position(weaken_pos);
			putative_island_length=0; //we do not have a putative_island any more
		}
	}
}

inline
double long KullbakCounter::EntropyThreshold
		(const vector<double> & bg_probabilities)
{
	double long entropy_threshold=0;
	unsigned int letters_in_alphabet=bg_probabilities.size()-1;
	double shifted_bg_probabilities[letters_in_alphabet];
	for (unsigned short removed_background_letter=1;
			removed_background_letter<=letters_in_alphabet;
			removed_background_letter++)
	{
		for (unsigned short j=1;j<=letters_in_alphabet;j++)
			shifted_bg_probabilities[j-1]=
					j==removed_background_letter
						?0
						:
						(bg_probabilities[j]/
							(1.-bg_probabilities[removed_background_letter]));
		double shifted_bg_entropy=0.;
		for (unsigned short j=1;j<=letters_in_alphabet;j++)
			if (j!=removed_background_letter)
				shifted_bg_entropy-=
						shifted_bg_probabilities[j-1]*log_2(shifted_bg_probabilities[j-1]);
		if (entropy_threshold<shifted_bg_entropy)
			entropy_threshold=shifted_bg_entropy;
	}
	//now, entropy_threshold is the maximum of all the shifted background
	//entropies
	return entropy_threshold;
}

inline
double long
KullbakCounter::EntropyDistanceFrom(const SymbolsCounter::PWM &pwm) const
{
	if (Symbols.pattern_length!=pwm.pattern_length)
		throw (* new SymbolsCounter::OtherLengthPWMException
				("Trying to measure distance btw a PWN and a KullbakCounter of different lengthes.\n"));
	//we suppose that both of them are created equal :)
	//I mean geometry
	double long Dist=0;
	for (unsigned int i=0;i<Symbols.pattern_length;i++)
	{
		if (Symbols.is_spaced && i==Symbols.spacer_5_end) i=Symbols.spacer_3_end+1;
		if (!if_strong_position[i]) continue;
		for (unsigned short j=1;j<=Symbols.letters_in_alphabet;j++)
			Dist+=Symbols.foreground_probability(i,j)*
							log_2(Symbols.foreground_probability(i,j)/pwm.foreground_probability(i,j));
	}
	return Dist;
}


#endif //_KULLBAK_COUNTER_HPP
