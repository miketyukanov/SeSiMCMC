/****************************************************************************\
SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2021
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For general description of the classes declared in the header, see headers.txt
$Id$
\****************************************************************************/

#include "ResultsSet.hpp"

char AtgcLetterOrder::letter(unsigned int index) const
{
	if (index==1) {return 'a';};
	if (index==2) {return 't';};
	if (index==3) {return 'g';};
	if (index==4) {return 'c';};
	return 'n';
}


unsigned int AtgcLetterOrder::atgcindex(unsigned int index)  const
{
	return index;
}

char AcgtLetterOrder::letter(unsigned int index)  const
{
	if (index==1) {return 'a';};
	if (index==2) {return 'c';};
	if (index==3) {return 'g';};
	if (index==4) {return 't';};
	return 'n';
}

unsigned int AcgtLetterOrder::atgcindex(unsigned int index)  const
{
	if (index==1) {return 1 /*'a'*/;};
	if (index==2) {return 4 /*'c'*/;};
	if (index==3) {return 3 /*'g'*/;};
	if (index==4) {return 2 /*'t'*/;};
	return 0;
}


Profile::Profile (const SequencesPile & sp, const MarkovChainState & mcs,
											const LookingForAtgcMotifsMultinomialGibbs & sampler,
											unsigned short was_length_adjusted_p,
											unsigned int minimal_motif_length,
											unsigned int maximal_motif_length,
											unsigned short was_mode_slow_p,
											retrieve_mode_type retrieve_mod,
											unsigned int retrieve_threshld,
											unsigned int trim_edges_par,
											unsigned int dumb_trim_edges_par,
											unsigned int common_background_par,
											const vector<double> & background_par,
											Diagnostics & diags,
											const LetterOrder & orderp,
											unsigned int flanks_len,
 											unsigned short show_lettrs_in_gap,
											const set <pair<unsigned int,double> > & G_values_p,
											unsigned int output_complements_as_in_seq):
												state(* new MarkovChainState(mcs)),
												counter(sampler.create_new_counter(sp,state)),
												letters_in_alphabet(4),
												min_motif_length(minimal_motif_length),
												max_motif_length(maximal_motif_length),
												spaced_sites(state.is_spaced()),
												spacer_5_end(state.spacer_5_end()),
												spacer_3_end(state.spacer_3_end()),
												sequences(sp.size()),
												was_length_adjusted(was_length_adjusted_p),
												was_mode_slow(was_mode_slow_p),
												min_seq_length(sp.min_length),
												max_seq_length(sp.max_length),
												retrieve_mode(retrieve_mod),
												retrieve_threshold(retrieve_threshld),
												flanks_length(flanks_len),
												show_letters_in_gap(show_lettrs_in_gap),
												trim_edges(trim_edges_par),
												dumb_trim_edges(dumb_trim_edges_par),
 												common_background(common_background_par),
												background(background_par),
												diagnostics(diags),
												order(orderp),
												G_values(G_values_p),
												output_complements_as_in_sequence
													(output_complements_as_in_seq),
												informative_part_length(state.informative_part_length()),
												shortest_significant_motif(sampler.shortest_significant_motif),
												SamplerType(typeid(sampler)),
												two_threads(
															SamplerType==typeid(LookingForAtgcMotifsInTwoThreadsMultinomialGibbs)||
															SamplerType==typeid(LookingForAtgcRepeatsInTwoThreadsMultinomialGibbs)
												),
												site_length(state.pattern_length()),
												site_length_before_trim(state.pattern_length())

{
	counter->calculate(state);
  //cout<< typeid(*counter).name() << endl;
	//	cout<<*counter<<flush;
	vector <double> scores(0);
	double threshold_score=0;
	vector <double> scores_by_seq(sp.size());
	sites_that_form_profile_count=counter->motifs();

	KullbakCounter * Kullbak=(KullbakCounter*)NULL;
	Kullbak=new KullbakCounter(*counter);

//KullbakCounter * Kullbak=(KullbakCounter*)(new KullbakCounter(*counter));
//	cout<<"\nSl="<<site_length<<"  fl="<<flanks_length<<"\n"<<flush;
	if
	(
		(trim_edges || dumb_trim_edges) &&
		SamplerType!=typeid(LookingForAtgcRepeatsInOneThreadMultinomialGibbs) &&
		SamplerType!=typeid(LookingForAtgcRepeatsInTwoThreadsMultinomialGibbs)
	)
	{
		unsigned int it_was_trimmed=0;

		if (trim_edges) //clever trimming
			Kullbak->weaken_outer_strong_islands(island_multiplier);

		if
		(
				SamplerType==typeid(LookingForAtgcPalindromesMultinomialGibbs) ||
				spaced_sites
		)
		//symmetric trimming
		{
			unsigned i=0;
			while
			(
					site_length>the_shortest_sensible_motif+i+i+state.spacer_length()
					&&
					! (Kullbak->if_a_strong_position[i])
					&&
					! (Kullbak->if_a_strong_position[site_length-1-i])
			)
			i++;
			if (i)
			{
				it_was_trimmed=1;
                try {
                    state.trim(i, i, sp);
                } catch (const DumbException &de) {
                    diagnostics << de.info;
                };
            }
		}
		else
		//asymmetric; spacer here is 0
		{
			unsigned i=0,j=0, something;
			do {
				something=0;
				if (
						site_length>the_shortest_sensible_motif+i+j
						&&
						! (Kullbak->if_a_strong_position[i])
				)
				{
					i++;
					something=1;
				}
				if (
						site_length>the_shortest_sensible_motif+i+j
						&&
						! (Kullbak->if_a_strong_position[site_length-1-j])
				)
				{
					j++;
					something=1;
				}
			} while (something);
			if (i || j)
			{
				it_was_trimmed=1;
                try {
                    state.trim(i, j, sp);
                } catch (const DumbException &de) {
                    diagnostics << de.info;
                };
            }
		}
		if (it_was_trimmed)
		{
			//redo the data initialization
			SymbolsCounter * counter_prev=counter;
			counter = sampler.create_new_counter(sp,state);
			counter->calculate(state);
			delete counter_prev;
			shortest_significant_motif=sampler.shortest_significant_motif;
			KullbakCounter * old_Kullbak= Kullbak;
			Kullbak= new KullbakCounter(*counter);
			delete old_Kullbak;
			site_length=state.pattern_length();
			informative_part_length=state.informative_part_length();
		}
	}

	vector <double> a(letters_in_alphabet*site_length);
	//it is matrix of q(pos,symbol)/p(symbol)
	//maybe later it will be split to a and p,
	//but now it is  not necessary
	//a(pos,symbol) is a[(symbol-1)*pattern_length+pos]
	//pos in pattern is 0-based, symbol is 1-based.
	for (unsigned short j=1;j<=letters_in_alphabet;j++)
		for (unsigned int i=0;i<site_length;i++)
			a[(j-1)*site_length+i]=
					(counter->foreground_probability(i,j)/
					 counter->background_probability(j)); //i is pos, j is symbol


	max_name_length=0;
	for (unsigned int p=0;p<sp.size();p++)
		if (sp.names[p].size()>max_name_length) max_name_length=
			sp.names[p].size();
//found the maximal name length
//
//

	for (unsigned int p=0;p<sp.size();p++)
	{
		string
			motif_text,
			leftFlank(flanks_length,' '),
			rightFlank(flanks_length,' ');

		if (!state.motif_present[p]) continue;

//		cout<<"   %"<<p<<"/"<<sp.size()<<";"<<sp[p].size()<<" "<<flush;
		unsigned int flank_5_border,flank_3_border;
		flank_5_border=
				(state.positions[p]>flanks_length)?
				state.positions[p]-flanks_length:0;

		flank_3_border=
				(state.positions[p]+site_length+flanks_length-1<sp[p].size()-1)?
				state.positions[p]+site_length+flanks_length-1:sp[p].size()-1;
//		cout<<" "<<flank_5_border<<":"<<state.positions[p]<<":"<<flank_3_border<<" "<<flush;
		if (two_threads	&& state.is_complement[p])
		{
			//we revert the complement contents to match the motif
			Atgc::complement(
						sp.nucleotides[p].substr(state.positions[p],site_length),
						motif_text
					);
//			cout<<"r"<<flush;
			for (unsigned int i=0;i<state.positions[p]-flank_5_border;i++)
				rightFlank[i]=
						Atgc::complement(sp.nucleotides[p][state.positions[p]-i-1]);
//			cout<<"5f"<<flush;
			for (unsigned int i=0;i<flank_3_border+1-site_length-state.positions[p];i++)
				leftFlank[flanks_length-i-1]=
						Atgc::complement(sp.nucleotides[p][state.positions[p]+site_length+i]);
//			cout<<"3f"<<flush;
		}
		else
		{
			motif_text=sp.nucleotides[p].substr(state.positions[p],site_length);
			for (unsigned int i=0;i<state.positions[p]-flank_5_border;i++)
				leftFlank[flanks_length-i-1]=
						sp.nucleotides[p][state.positions[p]-i-1];

//			cout<<"5f"<<flush;
			for (unsigned int i=0;i<flank_3_border+1-site_length-state.positions[p];i++)
				rightFlank[i]=
						sp.nucleotides[p][state.positions[p]+site_length+i];
//			cout<<"3f"<<flush;
		}
//		cout<<"*"<<flush;
		double score=0;
		vector<unsigned short>::const_iterator ptr=sp[p].begin()+state.positions[p];
		//ptr looks at first pattern letter
		for(unsigned int i=0;i<site_length;i++)
		{
			if (spaced_sites &&  i>=spacer_5_end && i<=spacer_3_end)
				continue;

      //for scoring, we revert the complement in any case.
			if (two_threads && state.is_complement[p])
			{
				unsigned short letter=Atgc::complement(ptr[site_length-1-i]);
				score+=log_2((double)a[(letter-1)*site_length+i]);
			}
			else
			{
				unsigned short letter=ptr[i];
				score+=log_2((double)a[(letter-1)*site_length+i]);
			}
		}

		for (unsigned int i=0;i<site_length;i++) //what position are to be capped
			motif_text[i]=Kullbak->if_a_strong_position[i]
										?
										toupper(motif_text[i])
										:
										tolower(motif_text[i]);

		Site aSite(sp.names[p],motif_text,leftFlank,rightFlank,state.positions[p],
				state.is_complement[p],score,0);
		/*cout<<"Adding Site("<<sp.names[p]<<","<<motif_text<<","<<leftFlank<<","<<rightFlank<<","<<state.positions[p]<<","<<
		state.is_complement[p]<<","<<score<<","<<0<<")"<<endl<<flush;
		if (find(aSite)!=end())
			cout<<"IT IS FOUND!!!"<<endl<<flush;
		cout<<"size before"<<size()<<" #"<<flush;
		*/
		insert(aSite);
		//cout<<"^ "<<"size after"<<size()<<endl<<flush;
		scores.push_back(score);
		scores_by_seq[p]=score;
	}


	if (retrieve_mode==yes_retrieve && retrieve_threshold==0)
		retrieve_mode=no_retrieve;

	if (retrieve_mode==yes_retrieve && retrieve_threshold>scores.size())
		retrieve_mode=no_retrieve;

	if (retrieve_mode==smart_retrieve && retrieve_threshold==0)
		retrieve_mode=pure_smart_retrieve;

	if (retrieve_mode==smart_retrieve && retrieve_threshold>scores.size())
		retrieve_mode=pure_smart_retrieve;

	if (retrieve_mode==yes_retrieve || retrieve_mode==smart_retrieve || retrieve_mode==pure_smart_retrieve)
	{
//	cout<<"-"<<flush;
		if (retrieve_mode!=pure_smart_retrieve)
		{
				sort(scores.begin(),scores.end());
				threshold_score=scores[retrieve_threshold-1];
		}
		//for pure_smart, we ise only local threshold
		for (unsigned int p=0;p<sp.size();p++)
		{
			if (site_length>sp[p].size()) continue;

			if (retrieve_mode==pure_smart_retrieve && (!state.motif_present[p])) continue;

			unsigned int leftmost_masked_position=0;
			for(unsigned int pos=0;pos<state.pattern_length();pos++)
				if(sp.mask[p][pos]&1u) leftmost_masked_position=pos;

			for (unsigned int pos=0;pos<sp[p].size()-site_length;pos++) //direct
			{

				if(sp.mask[p][pos+state.pattern_length()-1]&1u)
					leftmost_masked_position=pos+state.pattern_length()-1;
				if (leftmost_masked_position>=pos &&
						leftmost_masked_position<=pos+state.pattern_length()-1)
					continue;
				//if we are here, the tested position is unmasked
				double score=0;
				vector<unsigned short>::const_iterator ptr=sp[p].begin()+pos;
				//ptr looks at first pattern letter
				for(unsigned int i=0;i<site_length;i++)
				{
						if (spaced_sites && i==spacer_5_end) i=spacer_3_end+1;
						unsigned short letter=ptr[i];
						score+=log_2((double)a[(letter-1)*site_length+i]);
				}

				if (retrieve_mode!=pure_smart_retrieve && score<threshold_score) continue;
				//if we use yes_retrieve, it is the only test for this position
				//if we use yes_smart, we want the site to be stronger than
				//both its sequences site and the threshold,

				if
				(
					retrieve_mode!=yes_retrieve //i.e. one of smart modes
					&&
					state.motif_present[p]
					&&
					score<scores_by_seq[p]
				) continue;

				//DEBUG
				//if(retrieve_mode!=yes_retrieve )
				//{
				//	cout<<"DEBUG:"<<retrieve_mode<<"NAME="<<sp.names[p]<<"   score="<<score
				//		<<"   local thr="<<scores_by_seq[p]<<"   pres="<<state.motif_present[p]<<endl;
				//}
				//DEBUG
				string
					motif_text,
					leftFlank(flanks_length,' '),
					rightFlank(flanks_length,' ');

				motif_text=sp.nucleotides[p].substr(pos,site_length);

				unsigned int flank_5_border,flank_3_border;

				flank_5_border=
						(pos>flanks_length)?
						pos-flanks_length:0;

				flank_3_border=
						(pos+site_length+flanks_length-1<sp[p].size()-1)?
						pos+site_length+flanks_length-1:sp[p].size()-1;

				for (unsigned int i=0;i<pos-flank_5_border;i++)
					leftFlank[flanks_length-i-1]=
							sp.nucleotides[p][pos-i-1];

				for (unsigned int i=0;i<flank_3_border+1-site_length-pos;i++)
					rightFlank[i]=
							sp.nucleotides[p][pos+site_length+i];

				for (unsigned int i=0;i<site_length;i++) //what position are to be capped
					motif_text[i]=Kullbak->if_a_strong_position[i]
												?
												toupper(motif_text[i])
												:
												tolower(motif_text[i]);
				Site aSite(sp.names[p],motif_text,leftFlank,rightFlank,pos,0,score,1),
						 cSite(sp.names[p],motif_text,leftFlank,rightFlank,pos,0,score,0);
				if (find(cSite)==end())
					insert(aSite);
			}
			if (! two_threads) continue;

			//opposite strand

			leftmost_masked_position=0;
			for(unsigned int pos=0;pos<state.pattern_length();pos++)
				if(sp.mask[p][pos]&2u) leftmost_masked_position=pos;

			for (unsigned int pos=0;pos<sp[p].size()-site_length;pos++)
			{
				if(sp.mask[p][pos+state.pattern_length()-1]&2u)
					leftmost_masked_position=pos+state.pattern_length()-1;
				if (leftmost_masked_position>=pos &&
						leftmost_masked_position<=pos+state.pattern_length()-1)
					continue;
				//if we are here, the tested position is unmasked

				double score=0;
				vector<unsigned short>::const_iterator ptr=sp[p].begin()+pos;
				//ptr looks at first pattern letter
				for(unsigned int i=0;i<site_length;i++)
				{
					if (spaced_sites && i==spacer_5_end) i=spacer_3_end+1;
					unsigned short letter=Atgc::complement(ptr[site_length-1-i]);
					score+=log_2((double)a[(letter-1)*site_length+i]);
				}

				if (retrieve_mode!=pure_smart_retrieve && score<threshold_score) continue;
				//if we use yes_retrieve, it is the only test for this position
				//if we use yes_smart, we want the site to be stronger than
				//both its sequences site and the threshold,

				if
				(
					retrieve_mode!=yes_retrieve //i.e. one of smart modes
					&&
					state.motif_present[p]
					&&
					score<scores_by_seq[p]
				) continue;

				//it is a complenent, whatever

				string
					motif_text,
					leftFlank(flanks_length,' '),
					rightFlank(flanks_length,' ');

				Atgc::complement(
							sp.nucleotides[p].substr(pos,site_length),
							motif_text
						);

				unsigned int flank_5_border,flank_3_border;

				flank_5_border=
						(pos>flanks_length)?
						pos-flanks_length:0;

				flank_3_border=
						(pos+site_length+flanks_length-1<sp[p].size()-1)?
						pos+site_length+flanks_length-1:sp[p].size()-1;

				for (unsigned int i=0;i<pos-flank_5_border;i++)
					rightFlank[i]=
							Atgc::complement(sp.nucleotides[p][pos-i-1]);

				for (unsigned int i=0;i<flank_3_border+1-site_length-pos;i++)
					leftFlank[flanks_length-i-1]=
							Atgc::complement(sp.nucleotides[p][pos+site_length+i]);

				for (unsigned int i=0;i<site_length;i++) //what position are to be capped
					motif_text[i]=Kullbak->if_a_strong_position[i]
												?
												toupper(motif_text[i])
												:
												tolower(motif_text[i]);
				Site aSite(sp.names[p],motif_text,leftFlank,rightFlank,pos,1,score,1),
						 cSite(sp.names[p],motif_text,leftFlank,rightFlank,pos,1,score,0);
				if (find(cSite)==end())
					insert(aSite);
			}
		}
	}
	//DEBUG
	//cout<<"\n<P>@ "<<diags<<" @<P>\n";

}

ostream & Profile::text_output (ostream & o) const
{

	if (was_mode_slow)
	{
		o<<"#The best information content for different site lengthes:"<<endl;
		for
		(
			set <pair<unsigned int,double> >::iterator G_it=G_values.begin();
			G_it!=G_values.end();G_it++
		)
				o<<"# Lenght="<<G_it->first<<" G="<<G_it->second<<
				(G_it->first==site_length?" *** ":"")<<endl;
	}


	if (diagnostics.str()!="")
	{
		diagnostics.output_mode=comment_output;
			o<<diagnostics<<endl;
		diagnostics.output_mode=txt_output;
	}

	o<<"#The sampler was looking for";

	if (SamplerType==typeid(LookingForAtgcMotifsInTwoThreadsMultinomialGibbs))
		o<<" sites on both strands of ";
	if (SamplerType==typeid(LookingForAtgcMotifsInOneThreadMultinomialGibbs))
		o<<" sites on one strand of ";
	if (SamplerType==typeid(LookingForAtgcRepeatsInTwoThreadsMultinomialGibbs))
		o<<" direct repeat sites on both strands of ";
	if (SamplerType==typeid(LookingForAtgcRepeatsInOneThreadMultinomialGibbs))
		o<<" direct repeat sites on one strand of ";
	if (SamplerType==typeid(LookingForAtgcPalindromesMultinomialGibbs))
		o<<" palindomes in ";

	o<<sequences<<" sequences of length from "<<min_seq_length<<" through "
						<<max_seq_length<<"."<<endl;
	if (was_length_adjusted)
	{
		o<<"#the site length found is "<<site_length_before_trim<<endl;
		o<<"#the sampler has tested the lengths in ["<<min_motif_length<<","
				<<max_motif_length<<"]"<<endl;
	}
	else
	{
		o<<"#the site length found is "<<site_length_before_trim<<endl;
	}
	if (trim_edges || dumb_trim_edges)
	{
		o<<"#then the site was trimmed down to length "<<site_length<<endl;
	}
	if(site_length==min_motif_length && site_length != max_motif_length)
	{
		o<<"#!!! The motif length is the minimal tested length."<<endl;
		o<<"#!!! It possibly means that the result is almost garbage."<<endl;
		o<<"#!!! Try to shift the band of to use slow mode."<<endl;
	}
	if (spaced_sites)
		o<<"#the motif is spaced at ["<<spacer_5_end<<","
				<<spacer_3_end<<"]."<<endl;
	if (spaced_sites &&
			informative_part_length<=shortest_significant_motif)
	{
			o<<"#!!! The motif informative part is just the shortest"<<endl;
			o<<"#!!! significant pattern length."<<endl;
			o<<"#!!! It possibly means that the result is almost garbage."<<endl;
	}

	if (retrieve_mode!=no_retrieve)
	{
		o<<"#"<<"The motif is formed by "<<sites_that_form_profile_count<<" sites"<<endl;
		o<<"#"<<size()<<" motif sites altogether"<<endl;
		o<<"#"<<size()-sites_that_form_profile_count<<" were added by post-sampling processing."<<endl;
		o<<"#"<<"They are marked with \"*\" at line end."<<endl;
	}
	else
		o<<"#"<<size()<<" motif sites"<<endl;

	o<<"/Position are based on 0"<<endl;
/*
	if (output_additional_information)
	{
		o<<"<table><caption><i><font size=+1>The background symbols probabilities</font></i></caption>"<<endl;
		o<<"<tr><td width=75 >a</td><td width=75 >t</td><td width=75 >g</td><td width=75 >c</td></tr>"<<endl;
		for (unsigned int letter=1;letter<=letters_in_alphabet;letter++)
			o<<"<td>"<<setw(10)<<resetiosflags(ios::right)
				<<setiosflags(ios::left)<<setprecision(3)
				<<counter->background_probability(letter)<<"</td>";
		o<<"</tr></table>"<<endl;
		o<<"<br>"<<endl;
		o<<""<<"<br>"<<endl;
		o<<"<table>";
		o<<"<caption><i><font size=+1>The table of positional symbols probabilities forming the motif</font></i></caption>"<<endl;
		o<<"<tr><td width=75 >a</td><td width=75 >t</td><td width=75 >g</td><td width=75 >c</td></tr>"<<endl;
		for (unsigned int j=0;j<site_length;j++)
		{
			o<<"<tr>"<<endl;
			if (spaced_sites && j>=spacer_5_end && j<=spacer_3_end)
			{
				for (unsigned int letter=1;letter<=letters_in_alphabet;letter++)
					o<<"<td>"<<setw(10)<<resetiosflags(ios::right)
						<<setiosflags(ios::left)<<setprecision(3)
						<<counter->background_probability(letter)<<"</td>";
			}
			//o<<"<td>"<<setw(3)<<setiosflags(ios::right)<<j<<"</td>";
			else
			{
				for (unsigned int letter=1;letter<=letters_in_alphabet;letter++)
					o<<"<td>"<<setw(10)<<resetiosflags(ios::right)
						<<setiosflags(ios::left)<<setprecision(3)
						<<counter->foreground_probability(j,letter)<<"</td>";
			}
			o<<"</tr>"<<endl;
		}
		o<<"</table><hr>"<<endl;
	}

*/
	if (output_additional_information)
	{
		o<<"#The background model is ";
		if (!common_background)
			o<<"dynamic (it depends on positions of sites)";
		else if (!background.size())
			o<<"common, it is gathered from the data";
		else
			o<<"common, pre-given";
		o<<".\n";

		o<<"#The background symbols probabilities:"<<endl;
		//o<<"#  a         t         g         c"<<endl;
		o<<"#  "<<order.letter(1)<<"         "<<order.letter(2)<<"         "<<order.letter(3)<<"         "<<order.letter(4)<<endl;
		o<<"#  ";
		for (unsigned int letter=1;letter<=letters_in_alphabet;letter++)
		{
			o<<setw(10)<<resetiosflags(ios::right)
				<<setiosflags(ios::left)<<setprecision(3);
			/*if (common_background && background.size()==4)
					o<<background[order.atgcindex(letter)-1];
			else*/
					o<<counter->background_probability(order.atgcindex(letter));
		}
		o<<endl;
		o<<"#The positions symbols probabilities forming the motif:"<<endl;
		//o<<"#  a         t         g         c"<<endl;
		o<<"#  "<<order.letter(1)<<"         "<<order.letter(2)<<"         "<<order.letter(3)<<"         "<<order.letter(4)<<endl;
		for (unsigned int j=0;j<site_length;j++)
		{
			o<<"#  ";
			if (spaced_sites && j>=spacer_5_end && j<=spacer_3_end)
			{
				for (unsigned int letter=1;letter<=letters_in_alphabet;letter++)
				{
					o<<setw(10)<<resetiosflags(ios::right)
						<<setiosflags(ios::left)<<setprecision(3);
					/*if (common_background && background.size()==4)
						o<<background[order.atgcindex(letter)-1];
					else*/
						o<<counter->background_probability(order.atgcindex(letter));
				}
			}
			else
			{
				for (unsigned int letter=1;letter<=letters_in_alphabet;letter++)
					o<<setw(10)<<resetiosflags(ios::right)
						<<setiosflags(ios::left)<<setprecision(3)
						<<counter->foreground_probability(j,order.atgcindex(letter));
			}
			o<<endl;
		}
	}
	o<<resetiosflags(ios::right)<<resetiosflags(ios::left);

	o<<"#the output below is sites list. Every string is:"<<endl;
	o<<"#gene_id position length score ";
	if (two_threads)
		o<<"direction ";
	o<<"site ";
	if (two_threads&&
				output_complements_as_in_sequence) o<<"site-as-is";
	o<<endl;

	Profile::iterator siteit=begin();
	do
	{
		string id=siteit->GeneId;
		atgcstring Content_string=siteit->Content;

		if (spaced_sites)
			for(unsigned int i=spacer_5_end;i<=spacer_3_end;i++)
				Content_string[i]=SpacerSymbol;

		while (id.size()<max_name_length) id.append(" ");

		o<<id<<" \t";
		o<<setw(int(log10((float)max_seq_length))+3)<<siteit->position<<"\t";
		o<<setw(int(log10((float)max_seq_length))+3)<<site_length<<"\t";
		o<<setprecision(4)<<setw(10)<<siteit->score<<"\t";
		if (two_threads)
			o<<(siteit->is_complement ? "<":">");
		else o<<" ";
		o<<"\t";
		if (two_threads && output_complements_as_in_sequence)
		{
			if (siteit->is_complement)
				Content_string.out35(o);
			else
				o<<Content_string;
			o<<"\t";
		}
		o<<Content_string; //just a site output.
		if (retrieve_mode!=no_retrieve && siteit->is_retrieved_later) o<<"\t*";
		o<<endl;
		siteit++;
	} while (siteit!=end());
	o<<endl;
	return o;
}

ostream & Profile::html_output (ostream & o) const
{
	if (was_mode_slow)
	{
		o<<"The best information content for different site lengthes:";
		o<<"<br><br>"<<endl;
		o<<"<table border=1><tr><td>Length</td><td><G></td><td></td></tr>";
		for
		(
			set <pair<unsigned int,double> >::iterator G_it=G_values.begin();
			G_it!=G_values.end();G_it++
		)
				o<<"<tr><td>"<<G_it->first<<"</td><td>"<<G_it->second<<"</td>"
				<<(G_it->first==site_length?"<td>***</td>":"<td></td>")<<"</tr>"<<endl;
		o<<"<table><br>"<<endl;
	}
	if (diagnostics.str()!="")
	{
			o<<diagnostics<<endl;
	}

	o<<"The sampler was looking for ";

	if (SamplerType==typeid(LookingForAtgcMotifsInTwoThreadsMultinomialGibbs))
		o<<" sites on both strands of ";
	if (SamplerType==typeid(LookingForAtgcMotifsInOneThreadMultinomialGibbs))
		o<<" sites on one strand of ";
	if (SamplerType==typeid(LookingForAtgcRepeatsInTwoThreadsMultinomialGibbs))
		o<<" direct repeat sites on both strands of ";
	if (SamplerType==typeid(LookingForAtgcRepeatsInOneThreadMultinomialGibbs))
		o<<" direct repeat sites on one strand of ";
	if (SamplerType==typeid(LookingForAtgcPalindromesMultinomialGibbs))
		o<<" palindomes in ";

	o<<sequences<<" sequences of length from "<<min_seq_length<<" through "
						<<max_seq_length<<"."<<"<br><hr>"<<endl;

	if (was_length_adjusted)
	{
		o<<"The site length found is "<<site_length_before_trim<<".<br>"<<endl;
		o<<"The sampler has tested the lengths in ["<<min_motif_length<<","
				<<max_motif_length<<"]"<<".<br>"<<endl;
	}
	else
		o<<"The site length found is "<<site_length_before_trim<<".<br>"<<endl;

	if (trim_edges || dumb_trim_edges)
	{
		o<<"Then the site was trimmed down to length "<<site_length<<".<br>"<<endl;
	}

	if(site_length==min_motif_length && site_length != max_motif_length)
	{
		o<<"!!! The motif length is the minimal tested length."<<"<br>"<<endl;
		o<<"!!! It possibly means that the result is almost garbage."<<"<br>"<<endl;
		o<<"!!! Try to shift the band of to use slow mode."<<"<br>"<<endl;
	}
	if (spaced_sites)
		o<<"The motif is spaced at ["<<spacer_5_end<<","
				<<spacer_3_end<<"]."<<".<br>"<<endl;
	if (spaced_sites &&
			informative_part_length<=shortest_significant_motif)
	{
			o<<"!!! The motif informative part is just the shortest"<<"<br>"<<endl;
			o<<"!!! significant pattern length."<<"<br>"<<endl;
			o<<"!!! It possibly means that the result is almost garbage."<<"<br>"<<endl;
	}

	if (retrieve_mode != no_retrieve)
	{
		o<<"The motif is formed by "<<sites_that_form_profile_count<<" sites and ";
		o<<size()-sites_that_form_profile_count<<" additional sites found by postprocessing"<<"<br>"<<endl;
		o<<"The latter are marked with \"*\" in the table."<<"<br>"<<endl;
	}
	else
		o<<"The motif is formed by "<<size()<<" motif sites"<<".<br>"<<endl;
	o<<"<br>All the positions are 0-based.<br><hr>"<<endl;

	if (output_additional_information)
	{
		o<<"<p>The background model is ";
		if (!common_background)
			o<<"dynamic (it depends on positions of sites)";
		else if (!background.size())
			o<<"common, it is gathered from the data";
		else
			o<<"common, pre-given";
		o<<".</p>\n";
		o<<"<table><caption><i><font size=+1>The background symbols probabilities</font></i></caption>"<<endl;
		//o<<"<caption><i><font size=+1>The table of positional symbols probabilities forming the motif</font></i></caption>"<<endl;
		o<<"<tr><td width=75 >"<<order.letter(1)<<"</td><td width=75 >"<<order.letter(2)<<"</td><td width=75 >"<<order.letter(3)<<"</td><td width=75 >"<<order.letter(4)<<"</td></tr>"<<endl;
		for (unsigned int letter=1;letter<=letters_in_alphabet;letter++)
		{
			o<<"<td>"<<setw(10)<<resetiosflags(ios::right)
				<<setiosflags(ios::left)<<setprecision(3);
			/*if (common_background && background.size()==4)
				o<<background[order.atgcindex(letter)-1]<<"</td>";
			else*/
				o<<counter->background_probability(order.atgcindex(letter))<<"</td>";
		}
		o<<"</tr></table>"<<endl;
		o<<"<br>"<<endl;
		o<<""<<"<br>"<<endl;
		o<<"<table>";
		o<<"<caption><i><font size=+1>The table of positional symbols probabilities forming the motif</font></i></caption>"<<endl;
		//o<<"<tr><td width=75 >a</td><td width=75 >t</td><td width=75 >g</td><td width=75 >c</td></tr>"<<endl;
		o<<"<tr><td width=75 >"<<order.letter(1)<<"</td><td width=75 >"<<order.letter(2)<<"</td><td width=75 >"<<order.letter(3)<<"</td><td width=75 >"<<order.letter(4)<<"</td></tr>"<<endl;
		for (unsigned int j=0;j<site_length;j++)
		{
			o<<"<tr>"<<endl;
			if (spaced_sites && j>=spacer_5_end && j<=spacer_3_end)
			{
				for (unsigned int letter=1;letter<=letters_in_alphabet;letter++)
				{
					o<<"<td>"<<setw(10)<<resetiosflags(ios::right)
						<<setiosflags(ios::left)<<setprecision(3);
					/*if (common_background && background.size()==4)
						o<<background[order.atgcindex(letter)-1]<<"</td>";
					else*/
						o<<counter->background_probability(order.atgcindex(letter))<<"</td>";
				}
			}
			//o<<"<td>"<<setw(3)<<setiosflags(ios::right)<<j<<"</td>";
			else
			{
				for (unsigned int letter=1;letter<=letters_in_alphabet;letter++)
					o<<"<td>"<<setw(10)<<resetiosflags(ios::right)
						<<setiosflags(ios::left)<<setprecision(3)
						<<counter->foreground_probability(j,order.atgcindex(letter))<<"</td>";
			}
			o<<"</tr>"<<endl;
		}
		o<<"</table><hr>"<<endl;
	}
	o<<resetiosflags(ios::right)<<resetiosflags(ios::left);
	o<<"<table>"<<endl;
	o<<"<caption><i><font size=+1>The description of sites</font></i></caption>";
	o<<"<tr><td>Gene Id</td><td>Position</td><td>Length</td><td>Score</td>";
	if (two_threads)
		o<<"<td>Direction</td>";
	o<<"<td>Site in sequence</td>";
	if (two_threads && output_complements_as_in_sequence)
		o<<"<td width=10></td><td>Site-as-in-motif</td>";
	o<<"</tr>"<<endl;

	Profile::iterator siteit=begin();
	do
	{
		o<<"<tr>";
		string id=siteit->GeneId;

		o<<"<td>"<<id<<"</td>";
		o<<"<td align=center>"<<siteit->position<<"</td>";
		o<<"<td align=center>"<<site_length<<"</td>";
		o<<"<td align=center>"<<setprecision(4)<<siteit->score<<"</td>";
		if (two_threads)
			o<<"<td align=center>"<<(siteit->is_complement ? "&lt":"&gt")<<"</td>";
		if (two_threads&&output_complements_as_in_sequence)
		{
			o<<"<td><tt>";
			if (siteit->is_complement)
			{
				o<<"<font color=grey>";
				siteit->rightFlank.out35(o);
				o<<"</font>";
				for (unsigned int pos=siteit->Content.size();pos>0;pos--)
				{
					if (spaced_sites && (!show_letters_in_gap)
							&&  pos>=spacer_5_end+1 && pos<=spacer_3_end+1)
					{
						o<<SpacerSymbol;
						continue;
					}
					if (spaced_sites && show_letters_in_gap
							&&  pos==spacer_3_end+1)
						o<<"<font color=grey>";

					o<<Atgc::complement(siteit->Content[pos-1]);

					if (spaced_sites && show_letters_in_gap
							&&  pos==spacer_5_end+1)
						o<<"</font>";
				}
//		It was		siteit->Content.out35(o);
				o<<"<font color=grey>";
				siteit->leftFlank.out35(o);
				o<<"</font>";
			}
			else
			{
				o<<"<font color=grey>"<<siteit->leftFlank<<"</font>";
				for (unsigned int pos=0;pos<siteit->Content.size();pos++)
				{
					if (spaced_sites && (!show_letters_in_gap)
							&&  pos>=spacer_5_end && pos<=spacer_3_end)
					{
						o<<SpacerSymbol;
						continue;
					}
					if (spaced_sites && show_letters_in_gap
							&&  pos==spacer_5_end)
						o<<"<font color=grey>";

					o<<siteit->Content[pos];

					if (spaced_sites && show_letters_in_gap
							&&  pos==spacer_3_end)
						o<<"</font>";
				}
//It was	o<<siteit->Content();
				o<<"<font color=grey>"<<siteit->rightFlank<<"</font>";
			}
			o<<"</tt></td><td></td>";
		};
		o<<"<td><tt>"<<
				"<font color=grey>"<<siteit->leftFlank<<"</font>";

			for (unsigned int pos=0;pos<siteit->Content.size();pos++)
			{
				if (spaced_sites && (!show_letters_in_gap)
						&&  pos>=spacer_5_end && pos<=spacer_3_end)
				{
					o<<SpacerSymbol;
					continue;
				}
				if (spaced_sites && show_letters_in_gap
						&&  pos==spacer_5_end)
					o<<"<font color=grey>";

				o<<siteit->Content[pos];

				if (spaced_sites && show_letters_in_gap
						&&  pos==spacer_3_end)
					o<<"</font>";
			}
//It was	o<<siteit->Content();

		o<<"<font color=grey>"<<siteit->rightFlank<<"</font>"<<
				"</tt></td>"; //just site output
		if (retrieve_mode!=no_retrieve && siteit->is_retrieved_later)
			o<<"<td>*</td>";
		o<<"</tr>"<<endl;
		siteit++;
	} while (siteit!=end());
	o<<"</table><br>"<<endl;
	return o;
}

ostream & Profile::short_html_output (ostream & o) const
{
	o<<"<i><font size=+1>The list of sites</font></i>"<<"<br>"<<endl;
	o<<"<pre>"<<endl;
	Profile::iterator siteit=begin();
	do
	{
		string id=siteit->GeneId;
		atgcstring Content_string=siteit->Content;

		if (spaced_sites)
			for(unsigned int i=spacer_5_end;i<=spacer_3_end;i++)
				Content_string[i]=SpacerSymbol;

		while (id.size()<max_name_length) id.append(" ");

		o<<id<<"   "<<Content_string<<endl;
		siteit++;
	} while (siteit!=end());
	o<<"</pre>"<<endl;
	o<<"<br>"<<endl;
	return o;
}

ostream & Profile::xml_output (ostream & o, const string & id, const string & InputFileName) const
{
  o<<endl<<"  <motif length=\""<<site_length<<"\"";
	if (InputFileName!="") o<<" name=\""<<InputFileName<<"\"";
	if (id!="undefined") o<<" id=\""<<id<<"\"";
	o<<">"<<endl;

	o<<"    <word-list size=\""<<size()<<"\">\n";

	Profile::iterator siteit=begin();
	do
	{
		atgcstring Content_string=siteit->Content;

		/*if (spaced_sites)
			for(unsigned int i=spacer_5_end;i<=spacer_3_end;i++)
				Content_string[i]=XMLSpacerSymbol;
		*/
		o<<"      <word";
		o<<" name=\""<<siteit->GeneId<<"\"";
		o<<" score=\""<<setprecision(4)<<siteit->score<<"\"";
		o<<" length=\""<<site_length<<"\"";
		o<<" location=\""<<siteit->position+1<<"\"";
		o<<" strand=\""<<(siteit->is_complement ? "revcomp":"direct")<<"\">";
		o<<Content_string;
		o<<"</word>\n";
		siteit++;
	} while (siteit!=end());

	o<<"    </word-list>\n";


	o<<"  </motif>"<<endl;
/*  o<<"  <group";
	if (InputFileName!="") o<<" name=\""<<InputFileName<<"\"";
	o<<">"<<endl;
	o<<"  </group>"<<endl;
*/
	return o;
}



ostream & Profile::ibis_pwm_output (ostream & o, const string & TF, const string & suffix) const
{

		//o<<"#  a         t         g         c"<<endl;
		//o<<"#  "<<order.letter(1)<<"         "<<order.letter(2)<<"         "<<order.letter(3)<<"         "<<order.letter(4)<<endl;
    o<<">"<<TF<<" "<<TF<<suffix<<endl;	
    //o<<"#  "<<order.letter(1)<<"         "<<order.letter(2)<<"         "<<order.letter(3)<<"         "<<order.letter(4)<<endl;
    for (unsigned int j=0;j<site_length;j++)
    {
        if (spaced_sites && j>=spacer_5_end && j<=spacer_3_end)
        {
            for (unsigned int letter=1;letter<=letters_in_alphabet;letter++)
            {
                o<<setw(10)<<resetiosflags(ios::right)
                    <<setiosflags(ios::left)<<setprecision(3);
                /*if (common_background && background.size()==4)
                    o<<background[order.atgcindex(letter)-1];
                else*/
                    o<<counter->background_probability(order.atgcindex(letter));
            }
        }
        else
        {
            for (unsigned int letter=1;letter<=letters_in_alphabet;letter++)
                o<<setw(10)<<resetiosflags(ios::right)
                    <<setiosflags(ios::left)<<setprecision(3)
                    <<counter->foreground_probability(j,order.atgcindex(letter));
        }
        o<<endl;
    }
    o<<endl;
	o<<resetiosflags(ios::right)<<resetiosflags(ios::left);

	return o;
}

