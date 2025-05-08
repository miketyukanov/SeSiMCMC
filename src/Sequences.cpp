/****************************************************************************\
SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2021
$Id$
\****************************************************************************/

#include <stdio.h>
#include <string.h>
//mixed programming, alas....

#include <iostream>
#include <iomanip>
#include <string>

using namespace std;


#include "Sequences.hpp"
#include "Atgc.hpp"
#include "ResultsSet.hpp"

//#define __DEBUGLEVEL__ 5
//#define __DEBUGLEVEL__ 2

ostream & operator<< (ostream & os, const SequencesPile & sp)
{
	for (unsigned int i=0;i<sp.size();i++)
	{
	//write FastA record
		os<<">"<<sp.names[i]<<endl<<sp.nucleotides[i];
	//if (fseq.starred) os<<"*";//we use classic FastA with *
		os<<endl<<endl;
	}
	return os;
}

istream & operator>> (istream & is, SequencesPile & sp)
{
//it works with c-style reding. Let it be.
//Structure:
//Name (any string with > in first position)
//Sequence ..... (*-terminated or not)
//Name
//Sequence....
//

	string Buffer;

	string::size_type pos;

	enum  {skipbefore, name_read, name_and_NBRF_comment_read, reading_sequence,
									done} state=skipbefore;

	string sequence_text="", name="", next_name="";

	bool ANameIsReadForNextRecord=false;

	vector<unsigned short> seq;

	sp.Read_Diagnostics="";

	bool minus_mark_claimed=false;

	while (!is.eof() )
	{

		Buffer.erase();

		getline(is,Buffer);
		//getline(is,Buffer,'\n');

#if __DEBUGLEVEL__ >=5
		cerr<<"%"<<Buffer<<"%"<<"{"<<Buffer.size()<<"}"<<endl<<flush;
#endif
		if(Buffer.size()>0)
		{
			if (Buffer[Buffer.size()-1]=='\015')  //we procced possible trailing '^M'
				Buffer=Buffer.substr(0,Buffer.size()-1);
#if __DEBUGLEVEL__ >=5
			cerr<<"@"<<flush;
#endif
			if (Buffer[Buffer.size()-1]=='\012')  //we procced possible trailing '^J'
				Buffer=Buffer.substr(0,Buffer.size()-1);

#if __DEBUGLEVEL__ >=5
			cerr<<"@"<<flush;
#endif
			if (Buffer[Buffer.size()-1]=='\015')  //we procced possible trailing '^M'
				Buffer=Buffer.substr(0,Buffer.size()-1);
#if __DEBUGLEVEL__ >=5
			cerr<<"@"<<flush;
#endif
			if (Buffer[Buffer.size()-1]=='\012')  //we procced possible trailing '^J'
				Buffer=Buffer.substr(0,Buffer.size()-1);
		}	
#if __DEBUGLEVEL__ >=5
		cerr<<"@"<<endl<<flush;
#endif
#if __DEBUGLEVEL__ >=5
		cerr<<"*"<<Buffer<<"*"<<endl<<flush;
		cerr<<"#"<<state<<"#"<<endl<<flush;
#endif

		switch(state)
		{
		case skipbefore:
		case name_read:
			if (Buffer[0]=='>')
			{
				if (state==name_read)
				{
					sp.Read_Diagnostics+="A sequence name was overriden:\n";
					sp.Read_Diagnostics+=name;
					sp.Read_Diagnostics+="\nby:\n";
					sp.Read_Diagnostics+=Buffer.substr(1);
					sp.Read_Diagnostics+="\nmaybe, it is an input file inconsistency?\n";
				}
				name=Buffer.substr(1); //we omit the starting ">"
				state=name_read;
				break;

			};
			if (state==skipbefore && Buffer.find_first_not_of(" \t")!=string::npos)
                throw IOStreamException("Fasta reader has found nonblank line before name.\n");
            //it can be treat as error or as SequencesPile delimiter
        case name_and_NBRF_comment_read:
            //here, we have read the name from previous lines; the line can be
            //NBRF comment (it has obligatory "-") or
            //empty (omitted) or
            //or the sequence or its part.

            //the NBRF part can appear only once, so:
            if ((pos = Buffer.find_first_not_of(" \t")) != string::npos) {
                if (Buffer[pos] == '-') {
                    if (state==name_and_NBRF_comment_read)
                        throw IOStreamException(
                            "Fasta reader has found two NBRF comments in one record.\n");
                    state = name_and_NBRF_comment_read;
                    break; //from switch
                }
            } else
                break; //from switch

            //if we are here. the line is part of sequence or an error. So,
            //the sequence has been started

        case reading_sequence:
#if __DEBUGLEVEL__ >=5
			cerr<<"#reading-seq"<<flush;
#endif
			//now, we are reading a,t,g,c. First, let's reallocate fseq.sequence.
			if (state==reading_sequence) //we came here from switch
																	 //so, it makes sense to test whether
																	 //the line is empty or is it
																	 //the "> name" title line of next
																	 //FastA record in stream.
			{
				if (Buffer.find_first_not_of(" \t")==string::npos)
				{
					state=done;
					break;
				} //empty string. A sequence is over.
			}

			else state=reading_sequence; //we came here from previous case
																	 //because the line was not empty and not
																	 //NBRF comment

			if (Buffer[0]=='>')
			{
				next_name=Buffer.substr(1);
				ANameIsReadForNextRecord=true;
				state=done;
				break; //from switch
			};

			if (Buffer[Buffer.size()-1]=='*')  //we procced possible trailing '*'
															           //in FastA, it means end of record
			{
				state=done;
				Buffer=Buffer.substr(0,Buffer.size()-1);
			}

			if (is.eof()) state=done; //end of file, nothing more to do

			for (
						string::iterator Buf_it=Buffer.begin();
						Buf_it<Buffer.end();
						Buf_it++
					) while (isspace(*Buf_it)) Buffer.erase(Buf_it);

			sequence_text.append(Buffer);
			//iterator

#if __DEBUGLEVEL__ >=5
			cerr<<"#"<<endl<<flush;
#endif
		case done: //we never jump here :)
			break;

		}
		if (state==done)
		{
#if __DEBUGLEVEL__ >=5
			cerr<<"#done#"<<endl<<flush;
#endif
			//a record is read
			sp.names.push_back(name);
			try {
				sp.push_back(Atgc::string2atgc(sequence_text,seq));
            } catch (const AtgcException &A) {
                string message(A.info);
                message+="\nThe pre-processed string that gained the error was:\n";
				message+=sequence_text;
				message+="\n.\n";
                throw AtgcException(message.c_str());
            };
            sp.nucleotides.push_back(sequence_text);
            unsigned int current_length = sp.back().size();
            sp.caps.push_back(*new vector<unsigned short>(current_length));
			sp.mask.push_back(*new vector<unsigned short>(current_length));
			for (unsigned int pos=0;pos<current_length;pos++)
			{
				if( (tolower(sp.nucleotides.back()[pos])=='x')|| (tolower(sp.nucleotides.back()[pos])=='-'))
					sp.mask.back()[pos]=0x3;
				else
					sp.mask.back()[pos]=0;

				if( sp.nucleotides.back()[pos]=='-')
				{
					if (!minus_mark_claimed)
					{
						sp.Read_Diagnostics+="We treat \'-\' as a 1-base gap, which cannot be inside a site.\n";
						sp.Read_Diagnostics+="If you mean something else, change the FastA file.\n";
						minus_mark_claimed=true;
					}
				}
			}
			for (unsigned int pos=0;pos<current_length;pos++)
			{
				sp.caps.back()[pos]=0;
				if ( !sp.mask.back()[pos] && isupper(sp.nucleotides.back()[pos])
					) sp.caps.back()[pos]=1;
				//DeBUG-local
				//cerr<<"\n^^^"<<sp.names.back()<<" "<<pos<<" "<<sp.nucleotides.back()[pos]<<"  "<<sp.caps.back()[pos];
			};
#if __DEBUGLEVEL__ >=1
			cerr<<endl<<flush;
			cerr<<".................................................."<<endl<<flush;
			for (unsigned int kk=0;kk<current_length;kk++)
				cerr<<sp.back()[kk];
			cerr<<endl<<flush;
			for (unsigned int kk=0;kk<current_length;kk++)
				cerr<<sp.nucleotides.back()[kk];
			cerr<<endl<<flush;
			for (unsigned int kk=0;kk<current_length;kk++)
				cerr<<sp.mask.back()[kk];
			cerr<<endl<<flush;
			for (unsigned int kk=0;kk<current_length;kk++)
				cerr<<sp.caps.back()[kk];
			cerr<<endl<<".................................................."<<endl<<flush;
			cerr<<endl<<flush;
#endif
			if (sp.size()>1)
			{
				if (current_length>sp.max_length) sp.max_length=current_length;
				if (current_length<sp.min_length) sp.min_length=current_length;
			}
			else
				sp.min_length=sp.max_length=current_length;
			sequence_text="";
			state=skipbefore;
			if (ANameIsReadForNextRecord)
			{
				name=next_name;
				ANameIsReadForNextRecord=false;
				state=name_read;
			}
			//the last "if" works with the situation when there was no delimiter (empty line
			//or "*" after previous sequence.
			//So, the FastAInputStream object has stored the name started with ""
		}

#if __DEBUGLEVEL__ >=5
		cerr<<"*"<<Buffer<<"*"<<endl<<flush;
		cerr<<"$"<<state<<"$"<<endl<<flush;
#endif
	}
#if __DEBUGLEVEL__ >=5
	cerr<<"returning "<<endl<<flush;
#endif
	return is;
}

void SequencesPile::put_mask(const Profile & results,
		double masked_part)
//we rely of fact thet 0<masked_part<1
{
	unsigned int unmasked_wing,mask_size;
	mask_size=(int)((double)results.site_length*masked_part);

	if (mask_size==0) mask_size=1;
	if (mask_size%2 != results.site_length%2) mask_size--;
	//now, they are both odd or even together
	if (mask_size==0) mask_size=2;

	unmasked_wing=(results.site_length-mask_size)/2;

	Profile::iterator siteit=results.begin();

	do
	{
		unsigned int seq_no=0;
		while (names[seq_no]!=siteit->GeneId)
		{
			seq_no++;
			if (seq_no>=size())
			{
                throw DumbException("Try to put mask from results of an alien sequence.\n");
            };
        }

        unsigned int start_pos = siteit->position + unmasked_wing;
        unsigned int end_pos=
				siteit->position+results.site_length-unmasked_wing-1;
		if (!siteit->is_complement)
		{
			for (unsigned int pos=start_pos;pos<=end_pos;pos++)
				mask[seq_no][pos]|=0x1;
		}
		else
		{
			for (unsigned int pos=start_pos;pos<=end_pos;pos++)
				mask[seq_no][pos]|=0x2;
		}
		for (unsigned int pos=start_pos;pos<=end_pos;pos++)
			nucleotides[seq_no][pos]='x';
		siteit++;
	} while (siteit!=results.end());
};

