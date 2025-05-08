/****************************************************************************\
SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2021
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For general describtion of the classes declared in the header, see headers.txt
$Id$
\****************************************************************************/

#ifndef _MARKOV_CHAIN_STATE_HPP
#define _MARKOV_CHAIN_STATE_HPP

#include <string>
#include <vector>
#include <algorithm>

using namespace std;

#include "Sequences.hpp"
#include "Exception.hpp"
//this is for the_shortest_sensible_motif

const unsigned int the_shortest_sensible_motif=6;

class MarkovChainState
{
protected:
//all the global data

	unsigned int seqs;
	unsigned int length;

	unsigned int spaced;
	unsigned int spacer_5; //gap first position inside motif, 0-based
	unsigned int spacer_3; //gap last position motif, 0-based
	unsigned int spacer_len_asked; //it is the last gap length set by
	                         //SetSymmetricGap.
	//it can differ from spacer_3-spacer_5+1 that is real gap lenght.
	//   length            even   odd
	//  spacer_len_asked    (spacer_3-spacer_5+1)
	//      1                1     0
	//      2                1     2
	//      3                3     2
	//      4                3     4
	//      5                5     4
	//      6                5     6
	//      7                7     6

public:
//all the global data

	unsigned int sequences() const {return seqs;}   //const ref to seq
	unsigned int pattern_length() const {return length;}  //length

 	unsigned int is_spaced() const {return spaced;}   //spaced
	unsigned int spacer_5_end() const {return spacer_5;}   //spacer_5
	unsigned int spacer_3_end() const {return spacer_3;}   //spacer_3
	unsigned int spacer_length_asked() const {return spacer_len_asked;}
	//spacer_length aksed by last SetSymmetricGap
	void trim (unsigned int left, unsigned int right, const SequencesPile & sp)
	{
		if (spaced && (left != right))
		{
            throw DumbException("I cannot asymmetrically trim a spaced motif!\n");
        }
        if (length<the_shortest_sensible_motif+right+spacer_length())
		{
            throw DumbException("There was an attempt to trim down to a length\nthat is less than "
                                "the shortest sensible one!\n");
        }

        SetLength(length-left-right,sp,off); //caps_mode is not relevant here
		SetSymmetricGap(spacer_length());

		for (unsigned int i=0;i<seqs;i++)
		{
			if (!motif_present[i]) {continue;};
			if (!is_complement[i]) positions[i]+=left;
			else positions[i]+=right;
		};
	}

	unsigned int spacer_length() const
	{
		return spaced?(spacer_3-spacer_5+1):0;
	}

	unsigned int informative_part_length() const
	{
		return length-(spaced?(spacer_3-spacer_5+1):0);
	}
//seq-by-seq data
  vector<unsigned int> positions;
	vector<unsigned short> motif_present;
	vector<unsigned short> is_complement;
  vector<unsigned short> motif_was_too_long;
	//it is set to 1 if motif is too long
  //for its sequence or it is too long for each
  //caps sequence and the caps_mode is all

	MarkovChainState(unsigned int sequences_no):
			seqs(sequences_no),
			length(0),spaced(0),spacer_5(0),spacer_3(0),spacer_len_asked(0),
			positions(sequences_no),
			motif_present(sequences_no),is_complement(sequences_no),
			motif_was_too_long(sequences_no)
			{
				for (unsigned int i=0;i<seqs;i++)
				{
					positions[i]=0;//why not, and this position is valid for any seq
					motif_present[i]=1; //we suppose there is a motif by default;
															//Gibbs class can change it.
					is_complement[i]=0;
					motif_was_too_long[i]=0;
				};
			};

	void reset_positions()
	{
//it just resets all the positions, leaving the length_related things intact
		for (unsigned int i=0;i<seqs;i++)
		{
			positions[i]=0;//why not, and this position is valid for any seq
			if (!motif_was_too_long[i])
				motif_present[i]=1; //we suppose there is a motif by default;
													//Gibbs class can change it.
			is_complement[i]=0;
		};
	}

	void SetSymmetricGap(unsigned int gap)
	{
		spacer_len_asked=gap;
		if (gap>0)
		{
			//half_gap is the gap "wing" from middle to end
			if(length%2)  //length is even
			{
				unsigned int half_length=(length-1)/2,half_gap=(gap-1)/2;
				spacer_5=half_length-half_gap;
				spacer_3=half_length+half_gap;
			}
			else //length is odd
			{
				unsigned int half_length=length/2,half_gap=gap/2;
				spacer_5=half_length-half_gap;
				spacer_3=half_length+half_gap-1;
			}
			if (spacer_3>=spacer_5) spaced=1;
		}
		else
		{
			spaced=0;
			spacer_5=1;
			spacer_3=0;
		}
	};

	void SetLength(unsigned int len, const SequencesPile & sp, caps_mode caps)
	{
		length=len;
		SetSymmetricGap(spacer_len_asked);
		for (unsigned int p=0;p<sequences();p++)
		{
			unsigned int motif_does_fit=1;
			if (sp[p].size()<length) motif_does_fit=0;
			//it is too long for the sequence itself
			else
			{
				if (caps==all)
				//let's check whether there is at least one cap-island of necessary length
				{
					motif_does_fit=0;
					for (unsigned int pos=0;pos<sp[p].size()-length;pos++)
					{
						unsigned int all_caps=1;
						for (unsigned int int_pos=0;int_pos<length;int_pos++)
							if (! sp.caps[p][pos+int_pos] ) all_caps=0;
					  //if all_caps here =1, so we have an island
						if (all_caps)
						{
							motif_does_fit=1;
							break;
						}
					}
				}
			}
			if (! motif_does_fit )
			{
				motif_was_too_long[p]=1;
				motif_present[p]=0;
			}
			if (motif_was_too_long[p] && motif_does_fit)
			{
				motif_present[p]=0;
				//so it appears like absent and than can be
				//reanimated by LocalGibbsStep
				motif_was_too_long[p]=0;
			}
		}
	};
	//we set the pattern length and then reset the gap
	//using the last spacer_len_asked
	unsigned int useful_sequences() const
	{
		return seqs-count(motif_was_too_long.begin(),motif_was_too_long.end(),1);
	}

	unsigned int present_motifs() const
	{
		return useful_sequences()-count(motif_present.begin(),motif_present.end(),0);
	}

	vector<unsigned int>::reference operator[]
			(vector<unsigned int>::difference_type i)
	{
		return positions.operator[](i);
	}

	vector<unsigned int>::const_reference operator[]
			(vector<unsigned int>::difference_type i) const
	{
		return positions.operator[](i);
	}

	void Dump(ostream & os) const;

friend
	istream & operator>> (istream & in, MarkovChainState & mcs);
};

inline
unsigned int operator==(const MarkovChainState & mcs1,
		const MarkovChainState & mcs2)
{
	if (mcs1.sequences()!=mcs2.sequences()) return 0;
	if (mcs1.pattern_length()!=mcs2.pattern_length()) return 0;
	if (mcs1.is_spaced()!=mcs2.is_spaced()) return 0;
	if (mcs1.is_spaced())
	{
		if (mcs1.spacer_5_end()!=mcs2.spacer_5_end()) return 0;
		if (mcs1.spacer_3_end()!=mcs2.spacer_3_end()) return 0;
	}
	for (unsigned int i=0;i<mcs1.sequences();i++)
	{
		if(mcs1[i]!=mcs2[i]) return 0;
		if(mcs1.motif_present[i]!=mcs2.motif_present[i]) return 0;
		if(mcs1.is_complement[i]!=mcs2.is_complement[i]) return 0;
		if(mcs1.motif_was_too_long[i]!=mcs2.motif_was_too_long[i]) return 0;
	}
	return 1;
}

inline
unsigned int operator!=(const MarkovChainState & mcs1,
		const MarkovChainState & mcs2)
{
	return (mcs1==mcs2)?0:1;
}

//it is sclalar poroduct of two MCMC states.
//The old version (MarkovChainState.hpp <=v 1.9 ) was
//nucleotide-count-based StrictScalarProduct
//Now, we make it less strict.
inline
double operator*(const MarkovChainState & mcs1,
		const MarkovChainState & mcs2)
//it is a scalar product of two states.
//It counts motifs.
//Every existing motif add 1 to divider.
//If two motifs do not overlap,
//they add nothing to divident, while add 2 to divider.
//If their overlap is >= half of less motif, they add 2 to divider and to
//divident.
//If they overlap less, they add 2 to divider adn 1 to divident.
//We do not take gaps into account here.
{
	if (mcs2.sequences()!=mcs1.sequences()) return 0;
	unsigned int divident=0,divider=0;
	unsigned int radius=min(mcs1.pattern_length(),mcs2.pattern_length());
	for (unsigned int p=0;p<mcs1.sequences();p++)
	{
		if (!mcs1.motif_present[p]&&!mcs2.motif_present[p]) continue;
		//if we are here, at least one is present
		if (!mcs1.motif_present[p])
		{
			divider+=1;
			continue;
		}
		if (!mcs2.motif_present[p])
		{
			divider+=1;
			continue;
		}
		//here, the both are present
		if (mcs1.is_complement[p]!=mcs2.is_complement[p])
		{
			divider+=2;
			continue;
		}
		//here, the both are present and they have
		//we write here set-algorythms-based code
		//because it is the most general ans do not
		//suppose motif continuality,
		//and we suppose to use the operator not
		//too often, so we do not intend to optimise it
		//for speed.


		vector <unsigned int>
			positons_set1,
			positons_set2,
			positons_sets_intrs,
			positons_sets_union;

		for (unsigned int i=mcs1[p];i<mcs1[p]+mcs1.pattern_length();i++)
			positons_set1.push_back(i);

		for ( unsigned int i=mcs2[p];i<mcs2[p]+mcs2.pattern_length();i++)
			positons_set2.push_back(i);
		//of course, we can use more sophisticated shapes of pattern,
		//the only thing we need is to have sorted arrays, for
		//the using of "sets" is more expensive.
		set_intersection
			(
				positons_set1.begin(),
				positons_set1.end(),
				positons_set2.begin(),
				positons_set2.end(),
				back_inserter(positons_sets_intrs)
			);
		if (positons_sets_intrs.size()==0)
		{
			divider+=2;
			continue;
		}
		if (2*positons_sets_intrs.size()<radius)
		{
			divider+=2;
			divident+=1;
			continue;
		}
		divider+=2;
		divident+=2;
	}
	if (divider==0) return 0;
	else return double(divident)/double(divider);
}

inline
double StrictScalarProduct(const MarkovChainState & mcs1,
		const MarkovChainState & mcs2)
//it is a scalar product of two states.
//It counts bases.
//Every existing motif that do not overlap its cointerpart
//adds its length to divider.
//If two motifs do not overlap,
//they add nothing to divident.
//If they overlap, their uiuon is added to divider,
//and their intersection - to divident
{
	if (mcs2.sequences()!=mcs1.sequences()) return 0;
	double divident=0.,divider=0.;
	unsigned int pattern_lengths_sum=mcs1.pattern_length()+mcs2.pattern_length();
	for (unsigned int p=0;p<mcs1.sequences();p++)
	{
		if (!mcs1.motif_present[p]&&!mcs2.motif_present[p]) continue;
		//if we are here, at least one is present
		if (!mcs1.motif_present[p])
		{
			divider+=mcs2.pattern_length();
			continue;
		}
		if (!mcs2.motif_present[p])
		{
			divider+=mcs1.pattern_length();
			continue;
		}
		//here, the both are present
		if (mcs1.is_complement[p]!=mcs2.is_complement[p])
		{
			divider+=pattern_lengths_sum;
			continue;
		}
		//here, the both are present and they have
		//we write here set-algorythms-based code
		//because it is the most general ans do not
		//suppose motif continuality,
		//and we suppose to use the operatoe not
		//too often, so we do not intend to optimise it
		//for speed.

		vector <unsigned int>
			positons_set1,
			positons_set2,
			positons_sets_intrs,
			positons_sets_union;

		for (unsigned int i=0;i<mcs1.pattern_length();i++)
		{
			if (mcs1.is_spaced() && i==mcs1.spacer_5_end())
				i=mcs1.spacer_3_end()+1;
			positons_set1.push_back(mcs1[p]+i);
		}
		for (unsigned int i=0;i<mcs2.pattern_length();i++)
		{
			if (mcs2.is_spaced() && i==mcs2.spacer_5_end())
				i=mcs2.spacer_3_end()+1;
			positons_set2.push_back(mcs2[p]+i);
		}
		//of course, we can use more sophisticated shapes of pattern,
		//the only thing we need is to have sorted arrays, for
		//the using of "sets" is more expensive.
		set_union
			(
				positons_set1.begin(),
				positons_set1.end(),
				positons_set2.begin(),
				positons_set2.end(),
				back_inserter(positons_sets_union)
			);
		set_intersection
			(
				positons_set1.begin(),
				positons_set1.end(),
				positons_set2.begin(),
				positons_set2.end(),
				back_inserter(positons_sets_intrs)
			);
		divident+=positons_sets_intrs.size();
		divider+=positons_sets_union.size();
	}
	if (divider<=0) return 0;
	else return divident/divider;
}

inline
ostream & operator<< (ostream & o, const MarkovChainState & mcs)
{
	o<<"Motifs of length "<<mcs.pattern_length();
	if(mcs.is_spaced())
		o<<" with gap ["<<mcs.spacer_5_end()<<","<<mcs.spacer_3_end()<<"]";
	o<<" are located at:"<<endl;
	for (unsigned i=0;i<mcs.sequences();i++)
	{
		if (mcs.motif_present[i])
		{
			o<<mcs[i];
			if (mcs.is_complement[i])
				o<<"(compl)";
		}
		else
		{
			if (mcs.motif_was_too_long[i]>1)
				o<<"too-long";
			else
				o<<"absent";
		}
		o<<" ";
	}
	o<<endl;
	return o;
}

inline
void MarkovChainState::Dump(ostream & o) const
{
	o<<"Motifs of length "<<pattern_length();
	if(is_spaced())
		o<<" with gap ["<<spacer_5_end()<<","<<spacer_3_end()<<"]";
	o<<" are located at:"<<endl;
	for (unsigned i=0;i<sequences();i++)
	{
		if (motif_present[i])
		{
			o<<(*this)[i];
			if (is_complement[i])
				o<<"(compl)";
		}
		else
		{
			if (motif_was_too_long[i]>1)
				o<<"too-long("<<(*this)[i]<<")";
			else
				o<<"absent("<<(*this)[i]<<")";
		}
		o<<" ";
	}
	o<<endl;
}
inline
istream & operator>> (istream & in, MarkovChainState & mcs)
{
//We suppose that
//1)
//the first number on the first line of the
//record is the motif length,
//2)
//the first line starts with anything but the
//number, so the start can be delimiter (see 4)
//3)
//after it, the position information slots are given one-by-one
//delimited with any blank symbol.
//The normal position is a decimal number;
//a complemented is decimal number with text (no blanks!) containing
//triggraph "com"; the absent motif is marked with word "absent"
//4) Any space - delimited token differ from these three types
//means the end of MarkovChainState record. The end of line is
//omitted.
//
	string token;
	unsigned short bad_token;
	char * endp, *endp1;
	unsigned int k;
	//first, we clear it
	mcs.positions.clear();
	mcs.motif_present.clear();
	mcs.is_complement.clear();
	//now, we are reading tokens up to first unsigned int
	//which we believe to be pattern length.
	do
	{
		token="";
		in>>token;
		endp=NULL;
		k=strtoul(token.c_str(),&endp,10);
	} while( *endp && !cin.eof());
	if (token.empty())
	{
		//We did not find anything
        throw IOStreamException("Trying to read MarkovChainState from garbage.\n");
    } else {
        mcs.length=k;
		unsigned int g5,g3;
		do {
			char ch=in.get();
			if (ch=='\n') break;
			if (ch!='[') continue;
			token="";
			in>>token;
			endp=NULL;
			g5=strtoul(token.c_str(),&endp,10);
			if (endp && endp > token.c_str() && *endp==',')
			{
				g3=strtoul(endp+1,&endp1,10);
				if (endp1 && endp1 > endp+1 && *endp1==']' &&
						g5>0&&g3<mcs.pattern_length() && g5==mcs.pattern_length()-g3-1)
				{
					mcs.spacer_5=g5;
					mcs.spacer_3=g3;
					mcs.spaced=1;
					mcs.spacer_len_asked=g3-g5;
				}
			}
		} while(1);
    }
  //now, we start to read the positions, token-by-token
	mcs.seqs=0; //now, it is good tokens counter
	do
	{
		token="";
		bad_token=1;
		in>>token;

		if (token=="absent")
		{
			bad_token=0;
			mcs.seqs++;
			mcs.positions.push_back(0);
			mcs.is_complement.push_back(0);
			mcs.motif_present.push_back(0);
			mcs.motif_was_too_long.push_back(0);
			continue;
		};
		if (token=="too-long")
		{
			bad_token=0;
			mcs.seqs++;
			mcs.positions.push_back(0);
			mcs.is_complement.push_back(0);
			mcs.motif_present.push_back(0);
			mcs.motif_was_too_long.push_back(3);
			continue;
		};
		endp=NULL;
		k=strtoul(token.c_str(),&endp,10);
		if (*token.c_str() && !*endp) //just integer
		{
			bad_token=0;
			mcs.seqs++;
			mcs.positions.push_back(k);
			mcs.is_complement.push_back(0);
			mcs.motif_present.push_back(1);
			continue;
		}
		if
		(
				*token.c_str() && endp!=token.c_str() &&
				token.find("com",endp-token.c_str())!=string::npos
		)
		//integer with tail which include trigraph "com"
		{
			bad_token=0;
			mcs.seqs++;
			mcs.positions.push_back(k);
			mcs.is_complement.push_back(1);
			mcs.motif_present.push_back(1);
			continue;
		}
	} while (!bad_token && !in.eof());
	if (!mcs.sequences())
	{
		//We did not find anything
        throw IOStreamException("Empty MarkovChainState is suspisious.\n");
    }
    return in;
}
#endif //_MARKOV_CHAIN_STATE_HPP
