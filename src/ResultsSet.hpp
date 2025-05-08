/****************************************************************************\
SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2021
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For general description of the classes declared in the header, see headers.txt
$Id$
\****************************************************************************/

#ifndef _RESULTS_SET_HPP
#define _RESULTS_SET_HPP

#include <set>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <typeinfo>

using namespace std;

#include "Sequences.hpp"
#include "SymbolsCounter.hpp"
#include "KullbakCounter.hpp"
#include "MarkovChainState.hpp"
#include "MCMC.hpp"
#include "Diagnostics.hpp"

typedef enum {undef_retrieve=-1,no_retrieve=0,yes_retrieve=1,smart_retrieve=2,pure_smart_retrieve=3} retrieve_mode_type;

const char SpacerSymbol='.';
const char XMLSpacerSymbol='n';
const unsigned int island_multiplier=3;

struct LetterOrder
{
	virtual char letter (unsigned int index)  const = 0;
	virtual unsigned int atgcindex (unsigned int index) const = 0;
	//all indices are 1-based
	virtual ~LetterOrder() {};
};

struct AtgcLetterOrder : LetterOrder
{
	virtual char letter (unsigned int index) const;
	virtual unsigned int atgcindex (unsigned int index) const;
	//all indices are 1-based
};

struct AcgtLetterOrder : LetterOrder
{
	virtual char letter (unsigned int ATGCindex) const;
	virtual unsigned int atgcindex (unsigned int ATGCindex) const;
	//all indices are 1-based
};

inline
void output_text_header(ostream & out, const string & InputFileName)
{
	out<<"*Looking for DNA motifs Gibbs sampler (SeSiMCMC)"<<endl;
	out<<"#Input FastA file name:"<<endl;
	out<<"$"<<(InputFileName==""?"stdin":InputFileName)<<endl;
}

inline
void output_empty_text(ostream & out, const string & InputFileName)
{
	out<<"No reliable motifs were detected.\n";
}

inline
void output_html_header(ostream & out, const string & InputFileName)
{
	string IFN=(InputFileName==""?"stdin":InputFileName);
	out<<"<html>"<<endl
	   <<"<head>"<<endl
     <<"<meta http-equiv=\"Content-Type\" content=\"text/html\">"<<endl
     <<"<title>Yet Another Digging For DNA Motifs Gibbs Sampler"
		 <<" (SeSiMCMC) output.</title>"<<endl
		 <<"</head>"<<endl
     <<"<body bgcolor=\"white\">"
		 <<"<h2>Looking for DNA motifs Gibbs sampler (SeSiMCMC) result.</h2>"<<endl;
	out <<"<p>Input FastA file name: "<<IFN<<"</p>"<<endl;
	out<<"<!-- $"<<IFN<<" -->"<<endl;
}


inline
void output_empty_html(ostream & out)
{
	out<<"</h2>No reliable motifs were detected.</h2>"<<endl;
}

inline
void output_html_footer(ostream & out)
{
	out<<"</body></html>"<<endl;
}

inline
void output_xml_header(ostream & out)
{
	out<<"<?xml version=\'1.0\' encoding=\'UTF-8\'?>"<<endl
	<<"<smallbismark>"<<endl;
	out<<"<comment>All locations are 1-based</comment>"<<endl;
}

inline
void output_xml_footer(ostream & out)
{
	out<<"</smallbismark>"<<endl;
}



class flank:public string
{
public:
	flank(const string & s):string(s)
	{
		for (unsigned i=0;i<size();i++) at(i)=tolower(at(i));
	};
	ostream & out53(ostream &o) const;
	ostream & out35(ostream &o) const;
};

inline
ostream & operator << (ostream &o, const flank & fl)
{
	for (unsigned int i=0;i<fl.size();i++)
	{
		if (fl[i]!=' ')
			o<<fl[i];
		else
			o<<"&nbsp;";
	}
	return o;
}

inline
ostream & flank::out53(ostream &o) const
{
	for (unsigned int i=0;i<size();i++)
	{
		if (at(i)!=' ')
			o<<at(i);
		else
			o<<"&nbsp;";
	}
	return o;
}

inline
ostream & flank::out35(ostream &o) const
{
	for (unsigned int pos=size();pos>0;pos--)
	{
		char letter=at(pos-1);
		if (letter!=' ')
			o<<Atgc::complement(letter);
		else
			o<<"&nbsp;";
	}
	return o;
}

class atgcstring:public string
{
public:
	atgcstring(const string & s):string(s)
	{};
	ostream & out53(ostream &o) const;
	ostream & out35(ostream &o) const;
};

inline
ostream & atgcstring::out35(ostream &o) const
{
	unsigned int pos=0; //to awoid warnings about uninitialised use
	for (pos=size();pos>0;pos--)
	{
		char letter=0;
		letter=at(pos-1); //to awoid warnings about uninitialised use
		char letterc;
		try {
			letterc=Atgc::complement(letter);
        } catch (const AtgcException &a) {
            letterc = letter;
        }
        o<<letterc;
	}
	return o;
}

inline
ostream & atgcstring::out53(ostream &o) const
{
	for (unsigned int pos=0;pos<size();pos++)
		o<<at(pos);
	return o;
}

struct Site
{
	string GeneId;

	atgcstring Content;


	flank leftFlank, rightFlank;

	//content of a site is its "site content",
	//i.e the way that if form the motif ,
	//the flanks are compatible
	unsigned int position;
	unsigned short is_complement;
	double score;
	unsigned short is_retrieved_later;
	Site
		(
			const string & GI,
			const string & site,
			const string & lF,
			const string & rF,
			unsigned int pos,
			unsigned short is_c,
			double scor,
			unsigned short is_retrievd_later
		):GeneId(GI),Content(site),leftFlank(lF),rightFlank(rF),position(pos),is_complement(is_c),
			score(scor),is_retrieved_later(is_retrievd_later){};

};

inline
unsigned int operator== (const Site &S1, const Site & S2)
{
	if (S1.is_retrieved_later!=S2.is_retrieved_later) return 0;
	if (S1.GeneId!=S2.GeneId) return 0;
	if (S1.Content!=S2.Content) return 0;
//Actually, we can get equivaletnt Sites with different content string
//when they are created with different notation rules.....
	if (S1.leftFlank!=S2.leftFlank) return 0;
	if (S1.rightFlank!=S2.rightFlank) return 0;
	if (S1.position!=S2.position) return 0;
	if (S1.is_complement!=S2.is_complement) return 0;
	return 1;
};


inline
unsigned int operator< (const Site &S1, const Site & S2)
{
	if (S1.is_retrieved_later<S2.is_retrieved_later) return 1;
	if (S1.is_retrieved_later>S2.is_retrieved_later) return 0;
	if (S1.GeneId<S2.GeneId) return 1;
	if (S1.GeneId>S2.GeneId) return 0;
	if (S1.is_complement && !S2.is_complement) return 0;
	if (!S1.is_complement && S2.is_complement) return 1;
	if ((! S1.is_complement) && S1.position<S2.position) return 1;
	if ((! S1.is_complement) && S1.position>S2.position) return 0;
	if ((S1.is_complement) && S1.position<S2.position) return 0;
	if ((S1.is_complement) && S1.position>S2.position) return 1;
	if (S1.Content<S2.Content) return 1;
	if (S1.Content>S2.Content) return 0;
	if (S1.leftFlank<S2.leftFlank) return 1;
	if (S1.leftFlank>S2.leftFlank) return 0;
	if (S1.rightFlank<S2.rightFlank) return 1;
	if (S1.rightFlank>S2.rightFlank) return 0;
	return 0;
};

class Profile:public multiset< Site,less<Site> >
{
	//friend class Report;
	Profile();
	Profile(const Profile &);
	Profile & operator=(const Profile &);
	MarkovChainState & state;
	SymbolsCounter * counter;
	unsigned int sites_that_form_profile_count;
	unsigned int letters_in_alphabet;
	unsigned int min_motif_length;
	unsigned int max_motif_length;
	unsigned short spaced_sites;
	unsigned int spacer_5_end;
	unsigned int spacer_3_end;
	unsigned int sequences;
	unsigned short was_length_adjusted;
	unsigned short was_mode_slow;
	unsigned int min_seq_length;
	unsigned int max_seq_length;
	retrieve_mode_type retrieve_mode;
	unsigned int retrieve_threshold;
	unsigned int flanks_length;
	unsigned short show_letters_in_gap ; //default flanks
	unsigned int trim_edges;
	unsigned int dumb_trim_edges;
	unsigned int common_background;
	const vector<double> & background;
	Diagnostics & diagnostics;
	//equals to 0 if maximal chains to fail number is not overcomed,
	//otherwise equals to the number
	const LetterOrder & order;
	const set <pair<unsigned int,double> > & G_values;
	unsigned int output_complements_as_in_sequence;
	unsigned int informative_part_length;
	unsigned int shortest_significant_motif;
	//if 1, we output a 3'-5' complementary site as it is read from FastA,
	//not like it is accounted in motif. We output it as accounted in
	//next after "site-content" position of output,
	//which is intended for comments by April standard
	const type_info  & SamplerType;
	unsigned short two_threads;
	unsigned int max_name_length;
public:
	Profile (const SequencesPile & sp, const MarkovChainState & mcs,
		const LookingForAtgcMotifsMultinomialGibbs & sampler,
		unsigned short was_length_adjusted_p,
		unsigned int minimal_motif_length,
		unsigned int maximal_motif_length,
		unsigned short was_mode_slow_p,
		retrieve_mode_type retrieve_mode,
		unsigned int retrieve_threshold, //o if not ot retrieve
		unsigned int trim_edges_par,
		unsigned int dumb_trim_edges_par,
		unsigned int common_background_par,
		const vector<double> & background_par,
		Diagnostics & diags,
		const LetterOrder & orderp,
		unsigned int flanks_len=5,
		unsigned short show_lettrs_in_gap=0,
		const set <pair<unsigned int,double> > & G_values_p=
			*(new const set <pair<unsigned int,double> >),
		unsigned int output_complements_as_in_seq=1);

	unsigned int site_length;
	unsigned int site_length_before_trim;
	unsigned int output_additional_information;
	ostream & text_output (ostream & o) const;
	ostream & html_output (ostream & o) const;
	ostream & short_html_output (ostream & o) const;
	ostream & xml_output (ostream & o, const string & id, const string & InputFileName) const;
	ostream & ibis_pwm_output (ostream & o, const std::string & TF_name,const std::string & suffix) const;
};


#endif // _RESULTS_SET_HPP

