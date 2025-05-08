/****************************************************************************\
SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2021
$Id$
\****************************************************************************/

//#define __DEBUGLEVEL__ 5
//#define __DEBUGLEVEL__ 2

//here, we suppose that mask can be 2-bit:
//1 is direct strand masked; 2 is complement.
//Actually, now we have only 3 (masked) or 0
//(unmasked).

#include <math.h>

#include <iomanip>

#include <vector>

#include <iterator>

using namespace std;

#include "Exception.hpp"
#include "Sequences.hpp"
#include "MarkovChainState.hpp"
#include "SymbolsCounter.hpp"
#include "KullbakCounter.hpp"
#include "RandomMappings.hpp"
#include "MCMC.hpp"
#include "Atgc.hpp"

extern "C"{
	#include "Random.h"
}

const double pseudocounts_sum_on_annealing=1.5;

//we use mask[][] as test for masked positions

const char * core_git_id="$Id$";

//here, the random generator does exist and is initialised.

LookingForAtgcMotifsMultinomialGibbs::
	LookingForAtgcMotifsMultinomialGibbs
		(
			const SequencesPile & sp,
			double motif_absence_pr,
			unsigned long local_cycles,
			caps_mode initial_caps_opt,
			caps_mode annealing_caps_opt,
			caps_mode refinement_caps_opt,
			unsigned short if_adaptive_pseudocounts_opt,
			unsigned short if_common_background_opt,
			const vector<double> & background_opt,
			unsigned short be_quiet_opt,
			ostream & logstream
		):
	Y(2*(sp.max_length+1)),   //the size is enought even for compements counting
	Sequences(sp),
	letters_in_alphabet(4),
	sequences_count(sp.size()),
	current_sequence_index(0),
	be_quiet(be_quiet_opt),
	the_state (sp.size()),
	motif_absence_prior(motif_absence_pr),
	initial_caps(initial_caps_opt),
	annealing_caps(annealing_caps_opt),
	refinement_caps(refinement_caps_opt),
	if_adaptive_pseudocounts(if_adaptive_pseudocounts_opt),
	if_common_background(if_common_background_opt),
	background(background_opt),
	local_step_cycles_between_adjustments(local_cycles),
	overall_local_steps_made(0),
	overall_adjustments_made(0),
	log_stream(logstream),
	theState(the_state),
	permutation_forbidden(0)
{
#if __DEBUGLEVEL__ >=1
	be_quiet=0;
#endif
	//we cannot rely on minimal sequence length, for it can be
	//garbage like 2-based sequence.
	//
	//
	//default start pattern length
	//finding the shortest_significant_motif
	shortest_significant_motif=the_shortest_sensible_motif;

	default_pattern_lenght=
			min(shortest_significant_motif*6,Sequences.find_median_length()/5);

	default_pattern_lenght=max(default_pattern_lenght,shortest_significant_motif);

	//
	//longest_sensible_motif=(2*Sequences.min_length/3);


	//finding longest_sensible_motif
	longest_sensible_motif=(Sequences.find_median_length());


	if (motif_absence_prior>1) //we did not define it
		motif_absence_prior=0.1;


	if (local_step_cycles_between_adjustments==0)
	{
		local_step_cycles_between_adjustments=max((size_t)sp.max_length,sp.size());
		unsigned long lp=local_step_cycles_between_adjustments/10; //here is another parameter to play
		unsigned int order=1;
		while (lp>10) {lp=lp/10;order++;};
		if (lp <= 2 ) local_step_cycles_between_adjustments=2;
		else if (lp <= 5) local_step_cycles_between_adjustments=5;
		else local_step_cycles_between_adjustments=10; //"else" corresponds to
																							 //innermost "if", K&R
		for (;order>0;order--) local_step_cycles_between_adjustments=
															local_step_cycles_between_adjustments*10;
		//so, we have the integer from set 10,20,50,100,200,500,1000....
	}

	if(!be_quiet)
	{
		log_stream<<"We make "<<local_step_cycles_between_adjustments<<
		" local step cycles between positional adjustment"<<endl;
		log_stream<<"There are "<<sequences_count<<" sequences,";
		log_stream<<"The maximal sequence length is "<<Sequences.max_length<<endl;
		log_stream<<"The median sequence length is "<<Sequences.find_median_length()<<endl;
		log_stream<<"The minimal sequence length is "<<Sequences.min_length<<endl;
		log_stream<<"The default start motif length is "<<default_pattern_lenght<<endl;
		log_stream<<"The maximal reasonable motif length is "<<longest_sensible_motif<<endl;
		log_stream<<"The minimal reasonable motif length is "<<shortest_significant_motif<<endl;
		log_stream<<"Prior of motif absence	in a sequence is "<<motif_absence_prior<<" .\n";
		log_stream<<"Initial caps mode: ";
		switch (initial_caps)
		{
		case off:
			log_stream<<"off";break;
		case one:
			log_stream<<"each site contains at least one cap";break;
		case all:
			log_stream<<"each site contains caps only";break;
		default:log_stream<<"unknown";
		}
		log_stream<<".\n";

		log_stream<<"Annealing caps mode: ";
		switch (annealing_caps)
		{
		case off:
			log_stream<<"off";break;
		case one:
			log_stream<<"each site contains at least one cap";break;
		case all:
			log_stream<<"each site contains caps only";break;
		default:log_stream<<"unknown";
		}
		log_stream<<".\n";

		log_stream<<"Refinement caps mode: ";
		switch (refinement_caps)
		{
		case off:
			log_stream<<"off";break;
		case one:
			log_stream<<"each site contains at least one cap";break;
		case all:
			log_stream<<"each site contains caps only";break;
		default:log_stream<<"unknown";
		}
		log_stream<<".\n";

		//it is random in 0..max_0_based_pattern_position
	}
}

unsigned int LookingForAtgcMotifsMultinomialGibbs::reinit(unsigned int len)
{
	current_sequence_index=0;
	if(!be_quiet)
	{
		log_stream<<"We reinit the sampler with length "<<len<<" .\n";
	}
	return 0;
}

LookingForAtgcMotifsInOneThreadMultinomialGibbs::
	LookingForAtgcMotifsInOneThreadMultinomialGibbs
		(
			const SequencesPile & sp,
			double motif_absence_pr,
			unsigned long local_cycles,
			caps_mode initial_caps_opt,
			caps_mode annealing_caps_opt,
			caps_mode refinement_caps_opt,
			unsigned short if_adaptive_pseudocounts_opt,
			unsigned short if_common_background_opt,
			const vector<double> & background_opt,
			unsigned short be_quiet_opt,
			ostream & logstream
		):
	LookingForAtgcMotifsMultinomialGibbs
		(
			sp,
			motif_absence_pr,
			local_cycles,
			initial_caps_opt,
			annealing_caps_opt,
			refinement_caps_opt,
			if_adaptive_pseudocounts_opt,
			if_common_background_opt,
			background_opt,
			be_quiet_opt,
			logstream
		)
	{
		mode=one_thread;
		smode=no;
		if(!be_quiet) log_stream<<"We are looking for motifs on one thread.\n";
	};

LookingForAtgcMotifsInTwoThreadsMultinomialGibbs::
	LookingForAtgcMotifsInTwoThreadsMultinomialGibbs
		(
			const SequencesPile & sp,
			double motif_absence_pr,
			unsigned long local_cycles,
			caps_mode initial_caps_opt,
			caps_mode annealing_caps_opt,
			caps_mode refinement_caps_opt,
			unsigned short if_adaptive_pseudocounts_opt,
			unsigned short if_common_background_opt,
			const vector<double> & background_opt,
			unsigned short be_quiet_opt,
			ostream & logstream
		):
	LookingForAtgcMotifsMultinomialGibbs
		(
			sp,
			motif_absence_pr,
			local_cycles,
			initial_caps_opt,
			annealing_caps_opt,
			refinement_caps_opt,
			if_adaptive_pseudocounts_opt,
			if_common_background_opt,
			background_opt,
			be_quiet_opt,
			logstream
		)
	{
		mode=two_threads;
		smode=no;
		if(!be_quiet) log_stream<<"We are looking for motifs on both threads.\n";
	};

LookingForAtgcPalindromesMultinomialGibbs::
	LookingForAtgcPalindromesMultinomialGibbs
		(
			const SequencesPile & sp,
			double motif_absence_pr,
			unsigned long local_cycles,
			caps_mode initial_caps_opt,
			caps_mode annealing_caps_opt,
			caps_mode refinement_caps_opt,
			unsigned short if_adaptive_pseudocounts_opt,
			unsigned short if_common_background_opt,
			const vector<double> & background_opt,
			unsigned short be_quiet_opt,
			ostream & logstream
		):
	LookingForAtgcMotifsMultinomialGibbs
		(
			sp,
			motif_absence_pr,
			local_cycles,
			initial_caps_opt,
			annealing_caps_opt,
			refinement_caps_opt,
			if_adaptive_pseudocounts_opt,
			if_common_background_opt,
			background_opt,
			be_quiet_opt,
			logstream
		)
	{
		smode=palindromes;
		if(!be_quiet) log_stream<<"We are looking for palindromes.\n";
	};

LookingForAtgcRepeatsInOneThreadMultinomialGibbs::
	LookingForAtgcRepeatsInOneThreadMultinomialGibbs
		(
			const SequencesPile & sp,
			double motif_absence_pr,
			unsigned long local_cycles,
			caps_mode initial_caps_opt,
			caps_mode annealing_caps_opt,
			caps_mode refinement_caps_opt,
			unsigned short if_adaptive_pseudocounts_opt,
			unsigned short if_common_background_opt,
			const vector<double> & background_opt,
			unsigned short be_quiet_opt,
			ostream & logstream
		):
	LookingForAtgcMotifsInOneThreadMultinomialGibbs
		(
			sp,
			motif_absence_pr,
			local_cycles,
			initial_caps_opt,
			annealing_caps_opt,
			refinement_caps_opt,
			if_adaptive_pseudocounts_opt,
			if_common_background_opt,
			background_opt,
			be_quiet_opt,
			logstream
		)
	{
		mode=one_thread;
		smode=repeats;
		if(!be_quiet) log_stream<<"We are looking for double repeats on one thread.\n";
	};

LookingForAtgcRepeatsInTwoThreadsMultinomialGibbs::
	LookingForAtgcRepeatsInTwoThreadsMultinomialGibbs
		(
			const SequencesPile & sp,
			double motif_absence_pr,
			unsigned long local_cycles,
			caps_mode initial_caps_opt,
			caps_mode annealing_caps_opt,
			caps_mode refinement_caps_opt,
			unsigned short if_adaptive_pseudocounts_opt,
			unsigned short if_common_background_opt,
			const vector<double> & background_opt,
			unsigned short be_quiet_opt,
			ostream & logstream
		):
	LookingForAtgcMotifsInTwoThreadsMultinomialGibbs
		(
			sp,
			motif_absence_pr,
			local_cycles,
			initial_caps_opt,
			annealing_caps_opt,
			refinement_caps_opt,
			if_adaptive_pseudocounts_opt,
			if_common_background_opt,
			background_opt,
			be_quiet_opt,
			logstream
		)
	{
		mode=two_threads;
		smode=repeats;
		if(!be_quiet) log_stream<<"We are looking for double repeats on both threads.\n";
	};

SymbolsCounter* LookingForAtgcMotifsInOneThreadMultinomialGibbs::
	create_new_counter
	(
		const SequencesPile & sp,
		const MarkovChainState & mcs
	) const
{
	return new SymbolsCounter(sp,letters_in_alphabet,if_common_background,background);
};

SymbolsCounter* LookingForAtgcMotifsInTwoThreadsMultinomialGibbs::
	create_new_counter
	(
		const SequencesPile & sp,
		const MarkovChainState & mcs
	) const
{
	return new AtgcSymbolsCounter(sp,mcs,if_common_background,background);
}

SymbolsCounter* LookingForAtgcPalindromesMultinomialGibbs::
	create_new_counter
	(
		const SequencesPile & sp,
		const MarkovChainState & mcs
	) const
{
	return new SymmetricSymbolsCounter(sp,letters_in_alphabet,if_common_background,background);
};

SymbolsCounter* LookingForAtgcRepeatsInOneThreadMultinomialGibbs::
	create_new_counter
	(
		const SequencesPile & sp,
		const MarkovChainState & mcs
	) const
{
	return new DoubletSymbolsCounter(sp,letters_in_alphabet,if_common_background,background);
};

SymbolsCounter* LookingForAtgcRepeatsInTwoThreadsMultinomialGibbs::
	create_new_counter
	(
		const SequencesPile & sp,
		const MarkovChainState & mcs
	) const
{
	return new DoubletAtgcSymbolsCounter(sp,mcs,if_common_background,background);
}

double long LookingForAtgcMotifsMultinomialGibbs::find_maximum
		(
			unsigned int given_pattern_length,
			unsigned int minimal_pattern_length,
			unsigned int maximal_pattern_length,
			unsigned short adjust_pattern_length,
			unsigned short spaced,
			unsigned long steps_number_maximum_is_to_be_global_for,
			unsigned int annealings_number_maximum_is_to_be_global_for,
			unsigned int cycles_with_minor_change_to_say_cold,
			double minor_change_level,
			change_test_mode ctmode,
			unsigned int adjustments_during_annealing,
			unsigned int chain_fails_after,
			unsigned int chains_to_try,
			unsigned long cycles_per_annealing_attempt_par,
			unsigned int annealing_attempts,
			unsigned long time_limit,
			const LogRecorder & timer,
			//unsigned int annealings_with_different_length,
			//unsigned int reanneal_if_all_chains_fail,
			unsigned long max_steps
		) 
{

	if (time_limit && timer.how_long()>(double)time_limit)
        throw TimeLimitException(time_limit);

    if (max_steps == 0ul - 1) {
        max_steps = Sequences.size() * Sequences.find_median_length()
                    * Sequences.find_median_length() * Sequences.find_median_length();
    }

    if (local_step_cycles_between_adjustments < 5 * cycles_with_minor_change_to_say_cold)
        local_step_cycles_between_adjustments = 5 * cycles_with_minor_change_to_say_cold;
    //just to be sure we have enough time to anneal between
    //adjustments

    current_caps_mode = initial_caps;

    permutation_mode do_permutations = permut_adjustmentwise;

    if (permutation_forbidden)
        do_permutations = permut_never; //test workaround

    vector<unsigned int> schedule(Sequences.size());

    for (unsigned int i=0;i<Sequences.size();i++) schedule[i]=i;
	//initial shedule is how-it-was-given
	//ostream_iterator<unsigned int> out(cout," ");
	//copy(schedule.begin(),schedule.end(),out);
	//cout<<endl<<flush;

	unsigned long cycles_per_annealing_attempt;

	unsigned long max_annealing_steps=max_steps/2;
	if (!cycles_per_annealing_attempt_par)
	//we count the default value
		cycles_per_annealing_attempt=
				(unsigned long)(3.141275*Sequences.size()*
				Sequences.size()*Sequences.size());
	else
		cycles_per_annealing_attempt=cycles_per_annealing_attempt_par;

	unsigned int strong_positions_to_be_cold=the_shortest_sensible_motif-2;

	unsigned int minimal_length=minimal_pattern_length,
							 pattern_length=given_pattern_length,
							 maximal_length=maximal_pattern_length;
	SayLengthBounds(minimal_length,pattern_length, maximal_length);
	//here, we set the correct pattern length bounds
	//
//	pattern_length=pattern_length?pattern_length:default_pattern_lenght;

	//if pattern_length==0 we set the default
	//pattern_length=pattern_length<minimal_length?minimal_length:pattern_length;
	//pattern_length=pattern_length>maximal_length?maximal_length:pattern_length;
//it all now in SayLengthBounds

	unsigned long local_steps_made=0;
	unsigned long adjustments_made=0;
	unsigned long local_steps_made_after_adjustment=0;

	unsigned int initial_pattern_lenght=pattern_length;


	unsigned long steps_number_globalises_maximum=0;
	unsigned long annealing_steps_made=0;
	unsigned long local_cycles_with_minor_change=0;
	double motif_absence_prior_value=motif_absence_prior;
	unsigned short adjust_motif_length_on_annealing= 0;
	unsigned short adjust_motif_length_on_sampling= adjust_pattern_length;

	if(!be_quiet)
	{
		log_stream<<"\nAn annealing will go until "
			<<cycles_with_minor_change_to_say_cold
			<<" adjacent local cycles result in a minor change.\n";

		if (ctmode==site_positions)
			log_stream<<"A minor change is when	scalar products of states"<<endl
			<<"before and after the cycle is more than "
			<<minor_change_level<<".\n";
		if (ctmode==information)
			log_stream<<"A minor change is when	information distance between states"<<endl
			<<"before and after the cycle is less than "
			<<minor_change_level<<"\ntimes the information content of a state.\n";
		log_stream<<"Annealing attempt fails if it has not succeeded after "
			<<cycles_per_annealing_attempt<<" cycles."<<endl;
		log_stream<<"We will do "<<annealing_attempts<<" annealing attempts.\n";
		log_stream<<"The pattern length is"<<
				(adjust_motif_length_on_annealing?" ":" not")<<
				" changing during an annealing attempt.\n";
		log_stream<<"We look for motif between "<<minimal_length<<" and "<<
			maximal_length<<" length, starting with "<<pattern_length<<
			" .\n";
		if (adjust_pattern_length)
		{
			log_stream<<"We vary motif length during adjustments and "
					<<"from one attempt (chain) to another.\n";
		}
		else
		{
			log_stream<<"We do not vary motif length during adjustments and"
					<<"from one attempt to another.\n";
		}

		if (if_adaptive_pseudocounts)
		{
			log_stream<<"We reduce pseudocounts sum down to "<<pseudocounts_sum_on_annealing
					<<" during annealing.\n";
		}


	}

	//motif_absence_prior=0.;
	//we switch off theese facilities during annealing

	the_state.reset_positions();
	//later, the positiotions will be initialised correclty.
	//now, all we need is that all of ther are correct
	the_state.SetLength (pattern_length, Sequences,current_caps_mode);

	MarkovChainState state_before(sequences_count);
	state_before=theState;

	if (ctmode != site_positions && ctmode != information && ctmode != correlation)
	{
        throw DumbException("Unknown annealing test mode!!!.\n");
    }

    Symbols = create_new_counter(Sequences, theState);

    Symbols->calculate(theState);

    SymbolsCounter::PWM pwm_prev(*Symbols);
    //START_ANNEALING:
    double standard_B = Symbols->pseudocounts_sum();
    if (if_adaptive_pseudocounts)
        Symbols->change_pseudocounts_sum(pseudocounts_sum_on_annealing);

    if(!be_quiet) log_stream<<"Starting initial annealing attempts at "<<flush;

	unsigned int annealing_done=0;
	unsigned int annealing_attempt_no=0;
	unsigned long cycles_on_successful_attempt=0;
	while (!annealing_done) //annealing
	{
		unsigned long current_annealing_cycle=0;
		if (annealing_attempt_no)
		{
			annealing_attempt_no++;
			if (annealing_attempt_no>annealing_attempts)
			{
                throw AllInitialAnnealingAttemptsFailed(annealing_attempts, initial_pattern_lenght);
            }

            pattern_length = initial_pattern_lenght;

            if (!be_quiet)
                log_stream << "Initial annealing attempt #" << annealing_attempt_no << " at "
                           << flush;

            if (adjust_pattern_length) {
                //1 step uses initial_pattern_lenght
                //2 initial_pattern_lenght / sqrt(2)
                //3 initial_pattern_lenght * sqrt(2)
                //4 initial_pattern_lenght / 2
                //5 initial_pattern_lenght * 2
                if (annealing_attempt_no % 2) //3,5...
                {
                    for (unsigned int step_no = 5; step_no <= annealing_attempt_no; step_no += 4)
                        //5,9,...
								pattern_length*=2;
					if ((annealing_attempt_no%4)==3) //3,7,11
						pattern_length=(int)floor((double)pattern_length*sqrt((float)2));
					if (pattern_length>maximal_length) pattern_length=maximal_length;
                } else //2,4,6...
                {
                    double p_length = (double) pattern_length;
                    for (unsigned int step_no=4;
							step_no<=annealing_attempt_no;step_no+=4)
						//4,8,....
								p_length/=2;

					if ((annealing_attempt_no%4)==2) //2,6,10
						p_length/=sqrt((float)2);

					pattern_length=(int)floor(p_length);
					if (pattern_length<minimal_length) pattern_length=minimal_length;
                }
            }
        } else
            annealing_attempt_no++;

        annealing_steps_made = 0;

        the_state.SetLength (pattern_length, Sequences,current_caps_mode);
		if(!be_quiet) log_stream<<"length="<<pattern_length<<"...\n"<<flush;

		current_caps_mode=initial_caps;
		unsigned int sequences_long_enough=0;
		unsigned int sequences_OK=0;
		//initialization
		for (unsigned int i=0;i<sequences_count;i++)
		{
			//cout<<":"<<current_caps_mode<<" "<<initial_caps<<":"<<flush;
			if (theState.motif_was_too_long[i]) continue;
			sequences_long_enough++; //if we are here, the lenght is OK
			vector <unsigned int> allowed_positions;
			for (unsigned int pos=0;pos<Sequences[i].size()-the_state.pattern_length();pos++)
			{
				unsigned int pos_OK=0;
				if (current_caps_mode==off) pos_OK=1;
				else if (current_caps_mode==one)
				{
					for(unsigned int intenal_position=0;
						intenal_position<the_state.pattern_length();intenal_position++)
						if (Sequences.caps[i][pos+intenal_position]) pos_OK=1; //it is a cap here!!
				}
				else if (current_caps_mode==all)
				{
					pos_OK=1;
					for(unsigned int intenal_position=0;
						intenal_position<the_state.pattern_length();intenal_position++)
						if (! Sequences.caps[i][pos+intenal_position]) pos_OK=0; //nocap here!!
				}
				if (pos_OK) allowed_positions.push_back(pos);
			}
			//cout<<"AP["<<i<<"]:"<<allowed_positions.size()<<"   "<<flush;
			if (allowed_positions.size())
			{
				unsigned int mapped_position=
					(unsigned int)floor(uni()*(allowed_positions.size()));
				the_state[i]=allowed_positions[mapped_position];
				sequences_OK++;
				//cout<<the_state[i]<<"/"<<Sequences[i].size()<<endl<<flush;
			}
			else
				the_state.motif_present[i]=0;
		}
		if (sequences_long_enough*3<sequences_count)
            throw LenghtRequirementsFailed(sequences_OK, sequences_count);

        if (sequences_OK * 3 < sequences_count)
            throw TooResrtrictiveCapsMode(sequences_OK, sequences_count);

        if (!be_quiet)
            log_stream << "\nThe sampler initialised for the annealing attempt "
                       << annealing_attempt_no << endl
                       << theState << endl
                       << flush;

        Symbols = create_new_counter(Sequences, theState);
        //we are going to start sample the space with aStep() to anneal.
        //We can realise that the annealing is
        //over after cycles_with_minor_change_to_say_cold adjustments
        //have not change state too much.

        Symbols->calculate(theState);

        pwm_prev = *Symbols;

        current_caps_mode = annealing_caps;

        while (1) //annealing attempt
        {
            if (time_limit && timer.how_long() > (double) time_limit)
                throw TimeLimitException(time_limit);
            //if(!be_quiet) log_stream<<"\nTime check: limit="<<time_limit<<" timer.how_long()="<<timer.how_long()<<")...\n"<<flush;

            //		aStep();
            if (local_steps_made_after_adjustment
                < local_step_cycles_between_adjustments * sequences_count)
            //a local step
            {
                if (current_sequence_index == 0)
                //local steps cycles border, let's check the annealing-is-over
                //conditon
                {
                    current_annealing_cycle++;
                    if (current_annealing_cycle > cycles_per_annealing_attempt) {
                        if (!be_quiet) log_stream<<"Annealing attempt failed.\n"<<flush;
						break;
						//this annealing attempt failed
                    }
                    if (ctmode == site_positions)
					{
						double scal_prod=state_before*theState;
						if (scal_prod>=minor_change_level)
						{
							if ( ++local_cycles_with_minor_change>=
												 cycles_with_minor_change_to_say_cold)
							{
								annealing_done=1;
								cycles_on_successful_attempt=current_annealing_cycle;
								break;
							}
							//annealing is over, we say cold
						}
						else
							local_cycles_with_minor_change=0;
						state_before=theState;
					//if(!be_quiet) log_stream<<"Product="<<scal_prod<<endl<<flush;
					} //site_positions
					else //information
					{
						//cerr<<"PWM:  ";print_foreground_probablities(pwm_prev,cerr);
						//cerr<<"Symbols:  ";print_foreground_probablities(*Symbols,cerr);
						KullbakCounter Kullbak(*Symbols);
						double long distance = Kullbak.EntropyDistanceFrom(pwm_prev);
						if
							(
							 	Kullbak.strong_positions_counter>=strong_positions_to_be_cold &&
							 	Kullbak.strong_positions_counter>=pattern_length/4 &&
								distance/Kullbak.strong_positions_IC<=minor_change_level
							)
						{
							if ( ++local_cycles_with_minor_change>=
												 cycles_with_minor_change_to_say_cold)
							{
								annealing_done=1;
								cycles_on_successful_attempt=current_annealing_cycle;
								break;
							}
							//annealing is over, we say cold
						}
						else
							local_cycles_with_minor_change=0;
						pwm_prev=*Symbols;
						if(!be_quiet) log_stream<<"Distance(annealing)="<<distance<<
							"   StrongIC="<<Kullbak.strong_positions_IC<<"  in "<<
									Kullbak.strong_positions_counter<<" position."<<endl<<flush;
						//cerr<<"PWM after copy:  ";print_foreground_probablities(pwm_prev,cerr);
						//cerr<<"##################################################\n"<<flush;
					}
                }
                //if current sequence is too short,
                //we jump away from it.
                if (!theState.motif_was_too_long[schedule[current_sequence_index]]) {
                    unsigned short last_motif_protection=
							(theState.present_motifs()<1) ||
							(theState.present_motifs()==1 &&
							 theState.motif_present[schedule[current_sequence_index]]==1);
					double m_a_p_v=motif_absence_prior;
					if (last_motif_protection) motif_absence_prior=0;
	 		 		Symbols->exclude_sequence(theState,schedule[current_sequence_index]);
					LocalGibbsStep(schedule[current_sequence_index]);
					Symbols->include_sequence(theState,schedule[current_sequence_index]);
					if (last_motif_protection) motif_absence_prior=m_a_p_v;
					local_steps_made++;
					local_steps_made_after_adjustment++;
                }
                current_sequence_index++;
				current_sequence_index%=sequences_count;
				//current_sequence_index=(++current_sequence_index)%sequences_count;
			//set the next current sequence
            } else
            //an adjustment
            {
                if (do_permutations == permut_adjustmentwise)
                    RandomPermutationMapping(schedule);
                //				copy(schedule.begin(),schedule.end(),out);
                //				cout<<endl<<flush;

                if (adjustments_during_annealing) {
                    if (adjust_motif_length_on_annealing)
						PositionAndLengthAdjustment
								(minimal_length,maximal_length,
								 0, //no gap-test
								 local_steps_made,adjustments_made);
					//the previous line is just courtesy
					else PositionAdjustment(local_steps_made,adjustments_made,0 );
					//no gap-test
					adjustments_made++;
					local_steps_made_after_adjustment=0;
					NormaliseState();
					local_cycles_with_minor_change=0;
					pwm_prev=*Symbols;
                }
            }
//		we are here, so it was not a break "say cold";
//			if (++annealing_steps_made>max_annealing_steps) cout<<"RRRR"<<flush;
			if (++annealing_steps_made>max_annealing_steps) break;
        }
        // initial annealing is over if it was it was a break "say cold"
		// and so annealing_done=1;
	} //while(annealing_done)
	NormaliseState();
	if(!be_quiet) log_stream<<"\n\nThe sampler is initially annealed ... It took "<<cycles_on_successful_attempt<<" cycles on the attempt #"<<annealing_attempt_no<<"\n\n"<<flush;

	motif_absence_prior=motif_absence_prior_value;

	while(1)  //secondary annealing, still no length adjstms
	{
//		aStep();
if (time_limit && timer.how_long() > (double) time_limit)
    throw TimeLimitException(time_limit);

if (local_steps_made_after_adjustment < local_step_cycles_between_adjustments * sequences_count)
//a local step
{
    if (current_sequence_index == 0)
    //local steps cycles border, let's check the annealing-is-over
    //conditon
    {
        if (ctmode == site_positions) {
            double scal_prod = state_before * theState;
            if (scal_prod >= minor_change_level) {
                if (++local_cycles_with_minor_change >= cycles_with_minor_change_to_say_cold)
                    break;
                //annealing is over, we say cold
            } else
                local_cycles_with_minor_change = 0;
            state_before = theState;
            //if(!be_quiet) log_stream<<"Product="<<scal_prod<<endl<<flush;
        }
        //site_positions
        else //information
        {
            //cerr<<"PWM:  ";print_foreground_probablities(pwm_prev,cerr);
            //cerr<<"Symbols:  ";print_foreground_probablities(*Symbols,cerr);
            KullbakCounter Kullbak(*Symbols);
            double long distance = Kullbak.EntropyDistanceFrom(pwm_prev);
            if (Kullbak.strong_positions_counter >= strong_positions_to_be_cold
                && Kullbak.strong_positions_counter >= pattern_length / 4
                && distance / Kullbak.strong_positions_IC <= minor_change_level) {
                if (++local_cycles_with_minor_change >= cycles_with_minor_change_to_say_cold)
                    break;
                //annealing is over, we say cold
            } else
                local_cycles_with_minor_change = 0;
            pwm_prev = *Symbols;
            if (!be_quiet)
                log_stream << "Distance(annealing)=" << distance
                           << "   StrongIC=" << Kullbak.strong_positions_IC << "  in "
                           << Kullbak.strong_positions_counter << " position." << endl
                           << flush;
            //cerr<<"PWM after copy:  ";print_foreground_probablities(pwm_prev,cerr);
            //cerr<<"##################################################\n"<<flush;
        }
    }
    //if current sequence is too short,
    //we jump away from it.
    if (!theState.motif_was_too_long[schedule[current_sequence_index]]) {
        unsigned short last_motif_protection
            = (theState.present_motifs() < 1)
              || (theState.present_motifs() == 1
                  && theState.motif_present[schedule[current_sequence_index]] == 1);
        double m_a_p_v = motif_absence_prior;
        if (last_motif_protection)
            motif_absence_prior = 0;
        Symbols->exclude_sequence(theState, schedule[current_sequence_index]);
        LocalGibbsStep(schedule[current_sequence_index]);
        Symbols->include_sequence(theState, schedule[current_sequence_index]);
        if (last_motif_protection)
            motif_absence_prior = m_a_p_v;
        local_steps_made++;
        local_steps_made_after_adjustment++;
    }
    current_sequence_index++;
    current_sequence_index %= sequences_count;
    //current_sequence_index=(++current_sequence_index)%sequences_count;
    //set the next current sequence
} else
//an adjustment
{
    if (do_permutations == permut_adjustmentwise)
        RandomPermutationMapping(schedule);
    //copy(schedule.begin(),schedule.end(),out);
    //cout<<endl<<flush;

    if (adjustments_during_annealing) {
        if (adjust_motif_length_on_annealing)
            PositionAndLengthAdjustment(minimal_length,
                                        maximal_length,
                                        0, //no gap-test
                                        local_steps_made,
                                        adjustments_made);
        //the previous line is just courtesy
        else
            PositionAdjustment(local_steps_made, adjustments_made, 0);
        //no gap-test
        adjustments_made++;
        local_steps_made_after_adjustment = 0;
        NormaliseState();
        local_cycles_with_minor_change = 0;
        pwm_prev = *Symbols;
    }
}
//		aStep() is over;
		if (++annealing_steps_made>max_annealing_steps)
            throw LostInSpaceOnSecondaryAnnealing(max_annealing_steps, initial_pattern_lenght);
    }
    // secondary annealing is over

    if (if_adaptive_pseudocounts)
        Symbols->change_pseudocounts_sum(standard_B);
    NormaliseState();
	if(!be_quiet) log_stream<<"\n\nThe sampler is finally annealed ... \n\n"<<flush;

	//Now, we choose the more strong criteria for
	//steps_number_globalises_maximum :
	//from steps_number_maximum_is_to_be_global_for
	//or from annealings_number_maximum_is_to_be_global_for
	steps_number_globalises_maximum=
	(steps_number_maximum_is_to_be_global_for>
	 local_steps_made*annealings_number_maximum_is_to_be_global_for)?
	steps_number_maximum_is_to_be_global_for:
	local_steps_made*annealings_number_maximum_is_to_be_global_for;

	//Now, we want steps_number_maximum_is_to_be_global_for to be not
	//less than sequences_count*max_length
	//it is quite intuitive.
	//
	//Then (in the main cycle code), we will require to hold
	//the test "at least one adjustment was made after global
	//maximum" is better
	//

	steps_number_globalises_maximum=
	(steps_number_globalises_maximum>
	 	Sequences.max_length*sequences_count)?
	steps_number_globalises_maximum:
	 	Sequences.max_length*sequences_count;


	if(!be_quiet) log_stream<<"The maximum is global if "
			<<steps_number_globalises_maximum<<" steps do not change it.\n";


	unsigned long tracing_steps_made=0;

	unsigned long steps_after_maximum=0;

	MarkovChainState LastBestState(theState), TheBestStateAfterAdj(theState);

	double long F=0,maxG=0,G=0, bestLocalF=0;


	unsigned short an_adjustment_was_made_after_maximum=0;

	local_cycles_with_minor_change=0; //reset the annealing condition


	unsigned int adjacent_forbidden_adjustments=0;
	unsigned int failed_chains=0;
	double long maxGBeforeTheLastRegularAdjustment=0;
	unsigned short forbid_length_adjustment=0;
	unsigned long adjs_chain_fails_after=1000;
//we convert adjstmnt cycle into adjstmnt



	chains_to_try=chains_to_try?chains_to_try:3;
	//zero means auto-config, we choose 3


	if(!be_quiet) log_stream<<"A chain fails after "
			<<adjs_chain_fails_after<<" adjacent forbidden adjustments.\n";

	if(!be_quiet) log_stream<<"We have "
			<<chains_to_try<<" chain attempts.\n";

	MarkovChainState StateBeforeTheLastRegularAdjustment(theState);
  //we need it becouse an adj can get to strange result and all the chain
	//fails. We have point to rollback
	adjacent_forbidden_adjustments=0;
	if (do_permutations==permut_adjustmentwise) RandomPermutationMapping(schedule);
	//copy(schedule.begin(),schedule.end(),out);
	//cout<<endl<<flush;

	if (adjust_motif_length_on_sampling) G=PositionAndLengthAdjustment
						(minimal_length,maximal_length,spaced,
						 local_steps_made,adjustments_made);
	else
		G=PositionAdjustment(local_steps_made,adjustments_made,spaced);
//		G=PositionAdjustment(local_steps_made,adjustments_made,0);
	maxG=G;
	adjustments_made++;
	local_steps_made_after_adjustment=0;

	pwm_prev=*Symbols;

	current_caps_mode=refinement_caps;

	while(tracing_steps_made++<=max_steps )
	{
//		aStep();
if (time_limit && timer.how_long() > (double) time_limit)
    throw TimeLimitException(time_limit);

if (local_steps_made_after_adjustment < local_step_cycles_between_adjustments * sequences_count)
//a local step
{
    if (current_sequence_index == 0 && local_steps_made_after_adjustment >= sequences_count)
    //local steps cycles border, let's check the annealing-is-over
    //conditon (it is necessary for adjustment permission)
    {
        if (ctmode == site_positions) {
            double scal_prod = state_before * theState;
            if (scal_prod >= minor_change_level) {
                if (++local_cycles_with_minor_change >= cycles_with_minor_change_to_say_cold)
                    forbid_length_adjustment = 0;
                //annealing is over, we say cold
            } else {
                local_cycles_with_minor_change = 0;
                forbid_length_adjustment = 1;
            }
            state_before = theState;
            //if(!be_quiet) log_stream<<"Product="<<scal_prod<<endl<<flush;
        } //scalar
        else //informational
        {
            //cerr<<"PWM:  ";print_foreground_probablities(pwm_prev,cerr);
            //cerr<<"Symbols:  ";print_foreground_probablities(*Symbols,cerr);
            KullbakCounter Kullbak(*Symbols);
            double long distance = Kullbak.EntropyDistanceFrom(pwm_prev);
            if (Kullbak.strong_positions_counter >= strong_positions_to_be_cold
                && distance / Kullbak.strong_positions_IC <= minor_change_level) {
                if (++local_cycles_with_minor_change >= cycles_with_minor_change_to_say_cold)
                    forbid_length_adjustment = 0;
                //annealing is over, we say cold
            } else {
                local_cycles_with_minor_change = 0;
                forbid_length_adjustment = 1;
            }
            pwm_prev = *Symbols;

            if (!be_quiet)
                log_stream << "Distance(annealing)=" << distance
                           << "   StrongIC=" << Kullbak.strong_positions_IC << "  in "
                           << Kullbak.strong_positions_counter << " position." << endl
                           << flush;
            //cerr<<"PWM after copy:  ";print_foreground_probablities(pwm_prev,cerr);
            //cerr<<"##################################################\n"<<flush;
        }
    }
    //if current sequence is too short,
    //we jump away from it.
    if (!theState.motif_was_too_long[schedule[current_sequence_index]]) {
        Symbols->exclude_sequence(theState, schedule[current_sequence_index]);
        LocalGibbsStep(schedule[current_sequence_index]);
        Symbols->include_sequence(theState, schedule[current_sequence_index]);
        local_steps_made++;
        local_steps_made_after_adjustment++;
        if (bestLocalF < (F = Symbols->NegativeEntropy())) {
            bestLocalF = F;
            TheBestStateAfterAdj = theState;
            if (maxG
                < (G = Negative_Entropy_with_defined_patterns_per_pattern_position(theState, F))) {
                LastBestState = theState;
                maxG = G;
                if (!be_quiet)
                    log_stream << "***Local maximum (F=" << bestLocalF << ", G=" << maxG
                               << " ) : \n"
                               << theState << endl;
                steps_after_maximum = 0;
                an_adjustment_was_made_after_maximum = 0;
            }
            //whatever, this one is the best after the last adjsmt
        } else //the local maximum is still the best
        {
            if (++steps_after_maximum >= steps_number_globalises_maximum
                && an_adjustment_was_made_after_maximum) {
                //it is global!!!
                the_state = LastBestState;
                maxG = Negative_Entropy_with_defined_patterns_per_pattern_position(theState, F);
                NormaliseState();
                delete Symbols;
                if (!be_quiet)
                    log_stream << "***Global maximum, so returning ... (F=" << bestLocalF
                               << ", G=" << maxG << ")" << endl
                               << flush;
                return maxG;
            }
        }
    }
    current_sequence_index++;
    current_sequence_index %= sequences_count;
    //current_sequence_index=(++current_sequence_index)%sequences_count;
    //set the next current sequence
} else
//an adjustment
{
    if (do_permutations == permut_adjustmentwise)
        RandomPermutationMapping(schedule);
    //copy(schedule.begin(),schedule.end(),out);
    //cout<<endl<<flush;

    if (!forbid_length_adjustment)
    //adjustment permisson check
    {
        // optimisation, odnako... Posmotrim.....Nachinat Adjust s
        // nailuchsego a ne s tekushego.
        //			state_before=theState;
        //			the_state=TheBestStateAfterAdj;
        //			Symbols->calculate(theState);

        NormaliseState();
        StateBeforeTheLastRegularAdjustment = theState;
        maxGBeforeTheLastRegularAdjustment = maxG;
        adjacent_forbidden_adjustments = 0;
        if (adjust_motif_length_on_sampling)
            G = PositionAndLengthAdjustment(minimal_length,
                                            maximal_length,
                                            spaced,
                                            local_steps_made,
                                            adjustments_made);
        else
            G = PositionAdjustment(local_steps_made, adjustments_made, spaced);
        adjustments_made++;
        local_steps_made_after_adjustment = 0;
        an_adjustment_was_made_after_maximum = 1;
        //			Symbols->calculate(theState) was made by the adjustment;
        if (last_adjustment_has_changed_state) {
#if __DEBUGLEVEL__ >=1
					log_stream<<"Adjustment has changed state.\n"<<flush;
#endif
					bestLocalF=0;
					if(maxG<G)
					{
						LastBestState=theState;
						maxG=G;
						if(!be_quiet) log_stream<<"***Local maximum found by adjustment (G="<<maxG
								<<" )."<<endl;
						steps_after_maximum=0;
					}
					local_cycles_with_minor_change=0;
					pwm_prev=*Symbols;
        } else {
#if __DEBUGLEVEL__ >=1
					log_stream<<"Adjustment has not changed state.\n"<<flush;
#endif

				//  eto esli mi ego s nailuchshego nachinali.
				//	the_state=state_before;
				//	Symbols->calculate(theState);
				//
        }
#if __DEBUGLEVEL__ >=3
				log_stream<<"Dbg.\n"<<theState<<*Symbols<<endl<<flush;
#endif
    } else // adjustment is forbidden, the sampler is too hot.
    {
        adjacent_forbidden_adjustments++;
        if (!be_quiet)
            log_stream << "The sampler-is-annealed test failed (" << adjacent_forbidden_adjustments
                       << " 's adjacent one). \n";
        //       Commented forbiden adjustment!!!!!!
        //				if (adjustments_during_annealing)
        //					PositionAdjustment(local_steps_made,adjustments_made,spaced);
        //					PositionAdjustment(local_steps_made,adjustments_made,0);
        //       Commented forbiden adjustment!!!!!!
        else {
            if (!be_quiet)
                log_stream << "An adjustment is forbidden after " << local_steps_made
                           << " steps and " << adjustments_made << " regular adjustments made."
                           << "The state is :" << endl
                           << theState;
            if (!be_quiet)
                log_stream << "\n" << flush;
        }
        local_steps_made_after_adjustment = 0;
        local_cycles_with_minor_change = 0;
        pwm_prev = *Symbols;
        if (adjacent_forbidden_adjustments >= adjs_chain_fails_after) {
            if (!be_quiet)
                log_stream << "A chain has failed ";
            the_state = StateBeforeTheLastRegularAdjustment;
            Symbols->calculate(theState);
            maxG = maxGBeforeTheLastRegularAdjustment;
            an_adjustment_was_made_after_maximum = 0;
            steps_after_maximum = 0;
            failed_chains++;
            adjacent_forbidden_adjustments = 0;
            pwm_prev = *Symbols;
            local_cycles_with_minor_change = 0;
            forbid_length_adjustment = 0;
            if (!be_quiet)
                log_stream << "!!!\n";
            if (failed_chains >= chains_to_try) {
                if (!be_quiet)
                    log_stream << "All chains failed. Throwing exception at maxG=" << maxG << "!!!"
                               << endl
                               << flush;
                throw TooMuchChainsFailed(failed_chains, initial_pattern_lenght, maxG);
            }
        }
    }
}
//		aStep() is over;
	}
	if(!be_quiet) log_stream<<"lost in space. Throwing exception at maxG="<<maxG<<"!!!"<<endl<<flush;
    throw LostInSpaceOnTracing(max_steps, initial_pattern_lenght, maxG);
}

double LookingForAtgcMotifsMultinomialGibbs::find_maximum_slowly
		(
			unsigned int minimal_pattern_length,
			unsigned int maximal_pattern_length,
			unsigned short spaced,
			set <pair<unsigned int,double> > & G_values,
			unsigned long steps_number_maximum_is_to_be_global_for,
			unsigned int annealings_number_maximum_is_to_be_global_for,
			unsigned int cycles_with_minor_change_to_say_cold,
			double minor_change_level,
			change_test_mode ctmode,
			unsigned int adjustments_during_annealing,
			unsigned int chain_fails_after,
			unsigned int chains_to_try,
			unsigned long cycles_per_annealing_attempt,
			unsigned int annealing_attempts,
			unsigned long time_limit,
			const LogRecorder & timer,
			Diagnostics & diagnostics,
			unsigned long max_steps
		) 
{
	MarkovChainState MaxMCS(sequences_count);
	unsigned int at_least_one_succeeded=0;
	double long maximalG=0;

	unsigned int minimal_length=minimal_pattern_length,
							 maximal_length=maximal_pattern_length;

	unsigned int fake_pattern_length=0;

	SayLengthBounds(minimal_length,fake_pattern_length,maximal_length);

	float minor_change_test_ratio=1.;
	unsigned int globalise_ratio=1;
	unsigned int pattern_length=0;
	//cout << "min " <<minimal_length<<" max "<< maximal_length <<endl<<flush;
	for (
				pattern_length=minimal_length;
				pattern_length<=maximal_length;
				pattern_length++
			)
	{
		double long G=0.;
		unsigned int unreliable=0;
		//cout<<pattern_length<<endl<<flush;
		try
		{
			G=find_maximum
			(
				pattern_length,
				minimal_length, //do not make sense
			 	maximal_length, //do not make sense
				0,//adjust_pattern_length
				spaced,
				steps_number_maximum_is_to_be_global_for,
				annealings_number_maximum_is_to_be_global_for*globalise_ratio,
				cycles_with_minor_change_to_say_cold,
				minor_change_level/minor_change_test_ratio,
				ctmode,
				adjustments_during_annealing,
				chain_fails_after,
				chains_to_try,
				cycles_per_annealing_attempt,
				annealing_attempts,
				time_limit,
				timer
				//1 //reanneal_if_all_chains_failed
			);
			if(!be_quiet) log_stream<<"Ordinary return from find_maximum."<<endl<<the_state<<"G="<<G<<endl<<flush;
        } catch (const LookingForAtgcMotifsMultinomialGibbs::LostInSpaceOnSecondaryAnnealing &lo) {
            unreliable = 2;
            diagnostics<<lo;
			if(!be_quiet) log_stream<<"LostInSpaceOnSecondaryAnnealing caught."<<endl<<flush;
        } catch (const LookingForAtgcMotifsMultinomialGibbs::AllInitialAnnealingAttemptsFailed &lo) {
            unreliable = 2;
            diagnostics<<lo;
			if(!be_quiet) log_stream<<"AllInitialAnnealingAttemptsFailed caught."<<endl<<flush;
        } catch (const LookingForAtgcMotifsMultinomialGibbs::LostInSpaceOnTracing &lo) {
            unreliable = 1;
            diagnostics<<lo;
			G=lo.maxG;
			if(!be_quiet) log_stream<<"LostInSpaceOnTracing caught."<<endl<<the_state<<"G="<<G<<endl<<flush;
        } catch (const LookingForAtgcMotifsMultinomialGibbs::TooMuchChainsFailed &too) {
            unreliable = 1;
            diagnostics<<too;
			G=too.maxG;
			if(!be_quiet) log_stream<<"TooMuchChainsFailed caught."<<endl<<the_state<<"G="<<G<<endl<<flush;
        }
        //DEBUG

        if (!at_least_one_succeeded && unreliable < 2) //we started
		{
			if (!be_quiet) log_stream<<"&&& First length was counted. The best was:\n"<<theState<<"with G="<<G<<endl<<endl;
			MaxMCS=theState;
			maximalG=G;
			at_least_one_succeeded=1;
			G_values.insert(make_pair(pattern_length,G));
		}
		else
		{
			if (unreliable<2)
			{
				if (!be_quiet) log_stream<<"&&& A length was counted. The best was:\n"<<theState<<"with G="<<G<<endl<<endl;
				if (G>maximalG)
				{
					MaxMCS=theState;
					maximalG=G;
				};
				at_least_one_succeeded=1;
				G_values.insert(make_pair(pattern_length,G));
			}
			else
			unreliable=0;
		}
		//DEBUG
		//cout<<"<p>@<p>"<<diagnostics<<"<p>@<p>"<<flush;
	}
    if (!at_least_one_succeeded)
        throw SlowSearchFailed{};

    the_state.SetLength(MaxMCS.pattern_length(), Sequences, current_caps_mode);

    the_state = MaxMCS;

    Symbols=create_new_counter(Sequences,theState);
	//we are going to start sample the space with aStep() to anneal.
	//We can realise that the annealing is
	//over after cycles_with_minor_change_to_say_cold adjustments
	//have not change state too much.

	Symbols->calculate(theState);
	return maximalG;
}


void LookingForAtgcMotifsInOneThreadMultinomialGibbs::
		LocalGibbsStep(unsigned int sequence_no)
//it never calls Symbols.calculate(), for it is supposed to
//be called earlier or Symbols were not affected with
//something other than LocalGibbsStep.
{
	unsigned int outcome;
#if __DEBUGLEVEL__ >= 3
	double weights_sum, weights_sum_c, max_weight=0;
	unsigned int max_w_pos;
	if (theState.motif_was_too_long[sequence_no])
		log_stream<<"Trying to make local step for too short sequence "<<sequence_no<<endl;
#endif
	unsigned int motif_positions=
			Sequences[sequence_no].size()-the_state.pattern_length()+1;
	//the result we wait for is:
	//Y[0] is posterior for motif absence in sequence_no
	//Y[1..motif_positions] are posteriors for
	//motif presence in the position pos-1 (0-based)
	//
	//The local step is very similar for palindromes
	//and one thread.
	calculate_position_weights
			(
				the_state,
				sequence_no,
				motif_positions,
				*Symbols,
				Y
			);

#if __DEBUGLEVEL__ >= 3
	weights_sum=weights_sum_c=0;
	for (unsigned int i = 0 ;i<motif_positions;i++)
	{
		if (Y[i+1]>max_weight)
			max_weight=Y[(max_w_pos=i)+1];
		weights_sum+=Y[i+1];
	}
	log_stream<<"Local step, seq "<<
										sequence_no<<", "<<motif_positions
										<<" motif positions."<<endl;
	log_stream<<"max weight: Y["<<max_w_pos
			<<"]="<<max_weight<<endl<<"sum of all motifs:"<<weights_sum<<
			endl<<"motif absence weight is:"<<Y[0]<<endl;
#endif
	try{
		outcome=
				get_int_draw_from_a_ditribution
					(
						0,
						motif_positions,
						Y
					);
    } catch (const DumbException &de) {
        cerr << de << "\nExiting.\n";
        exit(1);
    }
#if __DEBUGLEVEL__ >= 3
	log_stream<<"Outcome is "<<outcome<<endl;
#endif
	if (outcome==0)
	{
		the_state.motif_present[sequence_no]=0;
		return;
	}
	//if (outcome<=motif_positions)
	{
		the_state.motif_present[sequence_no]=1;
		the_state[sequence_no]=outcome-1;
		return;
	}
}

void LookingForAtgcMotifsInTwoThreadsMultinomialGibbs::
		LocalGibbsStep(unsigned int sequence_no)
//it never calls Symbols.calculate(), for it is supposed to
//be called earlier or Symbols were not affected with
//something other than LocalGibbsStep.
{
	unsigned int outcome;
#if __DEBUGLEVEL__ >= 3
	double weights_sum, weights_sum_c, max_weight=0, max_weight_c=0;
	unsigned int max_w_pos,max_w_c_pos;
#endif
	unsigned int motif_positions=
			Sequences[sequence_no].size()-the_state.pattern_length()+1;
	//the result we wait for is:
	//Y[0] is posterior for motif absence in sequence_no read in direct way
	//Y[1..motif_positions] are posteriors for
	//motif presence in the position pos-1 (0-based) in thesame way
	//Y[motif_positions+1..2*motif_positions are posteriors of
	//motif present in position pos-motif_positions-1 (0-based) read as
	//complement
	//Y[2*motif_positions+1] is posterior of complementary sequence
	//and motif absence.
	calculate_position_weights
			(
				the_state,
				sequence_no,
				motif_positions,
				*Symbols,
				Y
			);

#if __DEBUGLEVEL__ >= 3
	weights_sum=weights_sum_c=0;
	for (unsigned int i = 0 ;i<motif_positions;i++)
	{
		if (Y[i+1]>max_weight)
			max_weight=Y[(max_w_pos=i)+1];
		weights_sum+=Y[i+1];
	}
	for (unsigned int i = 0 ;i<motif_positions;i++)
	{
		if (Y[motif_positions+1+i]>max_weight_c)
			max_weight_c=Y[(max_w_c_pos=i)+1+motif_positions];
		weights_sum_c+=Y[motif_positions+1+i];
	}
	log_stream<<"Local step with complement possibility, seq "<<
										sequence_no<<", "<<motif_positions<<" motif positions."<<endl;
	log_stream<<"Non-complemented max weight: Y["<<max_w_pos
			<<"]="<<max_weight<<endl<<"sum of all motifs:"<<weights_sum<<endl<<
			"motif absence weight is:"<<Y[0]<<endl;
	log_stream<<"Complemented max weight: Y["<<max_w_c_pos
			<<"]="<<max_weight_c<<endl<<"sum of all motifs:"<<weights_sum_c<<
			endl<<"motif absence weight is:"<<Y[2*motif_positions+1]<<endl;
#endif
	try{
		outcome=
				get_int_draw_from_a_ditribution
					(
						0,
						2*motif_positions+1,
						Y
					);
    } catch (const DumbException &de) {
        cerr << de << "\nExiting.\n";
        exit(1);
    }

#if __DEBUGLEVEL__ >= 3
	log_stream<<"Outcome is "<<outcome<<endl;
#endif
	unsigned int was_complement=the_state.is_complement[sequence_no];
	if (outcome==0)
	{
		the_state.motif_present[sequence_no]=0;
		the_state.is_complement[sequence_no]=0;
		if (was_complement!=0)
			(static_cast<AtgcSymbolsCounter*>(Symbols))->
					current_sequence_was_complemented(sequence_no);
		return;
	}
	if (outcome<=motif_positions)
	{
		the_state.motif_present[sequence_no]=1;
		the_state[sequence_no]=outcome-1;
		the_state.is_complement[sequence_no]=0;
		if (was_complement!=0)
			(static_cast<AtgcSymbolsCounter*>(Symbols))->
					current_sequence_was_complemented(sequence_no);
		return;
	}
	if (outcome<2*motif_positions+1)
	{
		the_state.motif_present[sequence_no]=1;
		the_state.is_complement[sequence_no]=1;
		the_state[sequence_no]=outcome-motif_positions-1;
		if (was_complement!=1)
			(static_cast<AtgcSymbolsCounter*>(Symbols))->
				current_sequence_was_complemented(sequence_no);
		return;
	}
	//if (outcome==2*motif_positions+1)
	{
		the_state.motif_present[sequence_no]=0;
		the_state.is_complement[sequence_no]=1;
		if (was_complement!=1)
			(static_cast<AtgcSymbolsCounter*>(Symbols))->
				current_sequence_was_complemented(sequence_no);
		return;
	}
}


void LookingForAtgcPalindromesMultinomialGibbs::
		LocalGibbsStep(unsigned int sequence_no)
//it never calls Symbols.calculate(), for it is supposed to
//be called earlier or Symbols were not affected with
//something other than LocalGibbsStep.
{
	unsigned int outcome;
#if __DEBUGLEVEL__ >= 3
	double weights_sum, weights_sum_c, max_weight=0;
	unsigned int max_w_pos;
#endif
	unsigned int motif_positions=
			Sequences[sequence_no].size()-the_state.pattern_length()+1;
	//the result we wait for is:
	//Y[0] is posterior for motif absence in sequence_no
	//Y[1..motif_positions] are posteriors for
	//motif presence in the position pos-1 (0-based)
	//
	//The local step is very similar for palindromes
	//and one thread.
	calculate_position_weights
			(
				the_state,
				sequence_no,
				motif_positions,
				*Symbols,
				Y
			);
#if __DEBUGLEVEL__ >= 3
	weights_sum=weights_sum_c=0;
	for (unsigned int i = 0 ;i<motif_positions;i++)
	{
		if (Y[i+1]>max_weight)
			max_weight=Y[(max_w_pos=i)+1];
		weights_sum+=Y[i+1];
	}
	log_stream<<"Local step, seq "<<
										sequence_no<<", "<<motif_positions
										<<" motif positions."<<endl;
	log_stream<<"max weight: Y["<<max_w_pos
			<<"]="<<max_weight<<endl<<"sum of all motifs:"<<weights_sum<<
			endl<<"motif absence weight is:"<<Y[0]<<endl;
#endif
	try{
		outcome=
				get_int_draw_from_a_ditribution
					(
						0,
						motif_positions,
						Y
					);
    } catch (const DumbException &de) {
        cerr << de << "\nExiting.\n";
        exit(1);
    }
#if __DEBUGLEVEL__ >= 3
	log_stream<<"Outcome is "<<outcome<<endl;
#endif
	if (outcome==0)
	{
		the_state.motif_present[sequence_no]=0;
		return;
	}
	//if (outcome<=motif_positions)
	{
		the_state.motif_present[sequence_no]=1;
		the_state[sequence_no]=outcome-1;
		return;
	}
}

//we count weights for all possible positions and gets one from them
//with probability proportional to the weight.


unsigned int LookingForAtgcMotifsMultinomialGibbs::
	ShiftOptimisation
	(
		MarkovChainState & state,
		unsigned int shift_magnitude, //supposed to be symmetric shift
		unsigned short & success,
		//success == 0 if all possibilities are masked
		long double & max_entropy,
		long double & G
	)
//it is an internal core procedure for all adjustment types
//we move state to shift_magnitude left and right,
//find the position with maximal NegativeEntropy
//and return the entropy in MaxNegativeEntropy, the optimal
//state in state and the return is boolean shows whether the
//adjustment has change state.
{
	double long F=0;
	unsigned int most_5_shift_abs=shift_magnitude,
							most_3_shift_abs=shift_magnitude;
	for (unsigned int i=0; i<sequences_count; i++)
	//seq-by-seq, including current
	{
		if (!state.motif_present[i])continue;
		if (!state.is_complement[i])
		{
			most_5_shift_abs=min(state[i],most_5_shift_abs);
			most_3_shift_abs=min((size_t)most_3_shift_abs,Sequences[i].size()-state[i]-
					state.pattern_length());
		}
		else //complement
		{
			most_3_shift_abs=min(state[i],most_3_shift_abs);
			most_5_shift_abs=min((size_t)most_5_shift_abs,Sequences[i].size()-state[i]-
					state.pattern_length());
		}
	}

	//	most_3_shift_abs and most_5_shift_abs are absolute modules of possible shifts

	MarkovChainState original_state (state);

#if __DEBUGLEVEL__ >= 2
	log_stream<<"Position optimization in [-"<<most_5_shift_abs<<".."
			<<most_3_shift_abs<<"] from state:\n"<<state;
#endif
	max_entropy=0;
	unsigned int max_entropy_shift=0;

	unsigned int last_found_masked_posintion_shift=Sequences.max_length+1;
	unsigned int last_found_uncapped_posintion_shift=Sequences.max_length+1;
	vector<unsigned int> last_found_capped_posintion_shifts(sequences_count);
	fill(last_found_capped_posintion_shifts.begin(),last_found_capped_posintion_shifts.end(),Sequences.max_length+1);

	for (unsigned int i=0; i<sequences_count; i++) //start
	//seq-by-seq, including current, we take up initial state
	{
		if (!state.motif_present[i]) continue;
		if (!state.is_complement[i])
			state[i]-=most_5_shift_abs;
		else
			state[i]+=most_5_shift_abs;
		for(unsigned int p=0;p<state.pattern_length()-1;p++)
		//we investigate the "primer" of length=pattern_length-1
		//then, we will start the "last-position-check" algorithm
		{
			if (!state.is_complement[i])
			{
				if (Sequences.mask[i][state[i]+p]&1u) last_found_masked_posintion_shift=p;
				if(current_caps_mode==all && Sequences.caps[i][state[i]+p]==0) last_found_uncapped_posintion_shift=p;
				if(current_caps_mode==one && Sequences.caps[i][state[i]+p]==1)
					last_found_capped_posintion_shifts[i]=p;
			}
			else
			{
				if(Sequences.mask[i][state[i]+state.pattern_length()-1-p]&2u) last_found_masked_posintion_shift=p;
				if(current_caps_mode==all && Sequences.caps[i][state[i]+state.pattern_length()-1-p]==0)
										last_found_uncapped_posintion_shift=p;
				if(current_caps_mode==one && Sequences.caps[i][state[i]+state.pattern_length()-1-p]==1)
					last_found_capped_posintion_shifts[i]=p;
			}
		}
	}

	unsigned short have_we_seen_an_allowed_position=0,
		was_previous_shift_calculated=0;

	success=0;


	//here, we just maximize the entropy factor F
	{
		SymbolsCounter * symb=create_new_counter(Sequences,state);
		for (unsigned int p=0; ; p++)
		{
			//we are now at (new) position p and we are testing the end of out mock motif,
			//i.e. the position p+length-1
			//for direct strand,
			//p is state[i];
			//new position is state[i]+length-1
			//step up is state[i}++
			//for complement,
			//the motif position is state[i]
			//the new position is state[i]
			//step up is state[i]--
			unsigned short the_current_is_forbidden=0;
			for (unsigned int i=0; i<sequences_count; i++) //start
			{
				if (!state.motif_present[i]) continue;
				if (!state.is_complement[i])
				{
					if(Sequences.mask[i][state[i]+state.pattern_length()-1]&1u)
					{
						last_found_masked_posintion_shift=p+state.pattern_length()-1;
					}
					if (current_caps_mode==all && Sequences.caps[i][state[i]+state.pattern_length()-1]==0)
					{
						last_found_masked_posintion_shift=p+state.pattern_length()-1;
					}
					if (current_caps_mode==one)
					{
						if(Sequences.caps[i][state[i]+state.pattern_length()-1]==1)
							last_found_capped_posintion_shifts[i]=p+state.pattern_length()-1;
					}
				}
				else
				{
					if(Sequences.mask[i][state[i]]&2u)
					{
						last_found_masked_posintion_shift=p+state.pattern_length()-1;
					}
					if (current_caps_mode==all && Sequences.caps[i][state[i]]==0)
					{
						last_found_uncapped_posintion_shift=p+state.pattern_length()-1;
					}
					if (current_caps_mode==one)
					{
						if(Sequences.caps[i][state[i]]==1)
							last_found_capped_posintion_shifts[i]=p+state.pattern_length()-1;
					}
				}
			}


			if (last_found_masked_posintion_shift>=p &&
					last_found_masked_posintion_shift<=p+state.pattern_length()-1)
						the_current_is_forbidden=1;
			else
			{

				if (current_caps_mode==all &&
					(
						last_found_uncapped_posintion_shift>=p
						&&
						last_found_uncapped_posintion_shift<=p+state.pattern_length()-1
					)
				)
				the_current_is_forbidden=1;
				//it is masked because of an older column.

				if (current_caps_mode==one)
				{
					for (unsigned int i=0; i<sequences_count; i++) //start
					{
						unsigned int this_seq_OK=0;
						if (
								last_found_capped_posintion_shifts[i]>=p
								&&
								last_found_capped_posintion_shifts[i]<=p+state.pattern_length()-1
							) //there is a cap inside
						{
							this_seq_OK=1;
						}
						if (!this_seq_OK) the_current_is_forbidden=1;
					}
				}
			}
			if (!the_current_is_forbidden)
			{
				if (!was_previous_shift_calculated)
				{
					symb->calculate(state);
					F=symb->NegativeEntropy();
				}
				else
				{
					symb->calculate_incrementally_after_3_end_shift(state);
					F+=symb->NegativeEntropy_change_after_3_end_shift();
#if __DEBUGLEVEL__ >= 1
					SymbolsCounter * ac=create_new_counter(Sequences,state);
					ac->calculate(state);
					if (!(*ac==*symb)) log_stream<<"Incremental 3'-shift/nonincremental discrepancy!!!\n(In PositionAdjustment)\n"<<
								"state is: /n"<<state<<
								"\ninc:\n"<<*symb<<"\ndirect:\n"<<*ac<<"\n"<<flush;
					double long ent=ac->NegativeEntropy();
					if (fabs(F-ent) > 1E-9) log_stream<<"Incremental entropy discrepancy:\n(In PositionAdjustment)\n"<<
						endl<<"increment is "<<F<<" direct is "<<ent<<endl<<flush;
					delete ac;
#endif
				}
				if (F>max_entropy || (!have_we_seen_an_allowed_position))
				{
					max_entropy=F;
					max_entropy_shift=p;
				}
				have_we_seen_an_allowed_position=1;
				was_previous_shift_calculated=1;
#if __DEBUGLEVEL__ >= 1
				log_stream<<(int)p-(int)most_5_shift_abs<<"  "<<F<<endl;
#endif
			}
			else
			{
				was_previous_shift_calculated=0;
#if __DEBUGLEVEL__ >= 1
				log_stream<<(int)p-(int)most_5_shift_abs<<" masked"<<endl;
#endif
			}
			if (p == most_5_shift_abs+most_3_shift_abs) break;
			//it was exit condition of the cycle
			for(unsigned int ic=0; ic<sequences_count; ic++) //3'-shift
			//seq-by-seq, including current
			{
				if (!state.motif_present[ic]) continue;
				if (!state.is_complement[ic])
					state[ic]++;
				else
					state[ic]--;
			}
		}
		delete symb;
	}
#if __DEBUGLEVEL__ >= 2
	log_stream<<"The best shift is "<<
			(int)max_entropy_shift-(int)most_5_shift_abs<<" with F="<<max_entropy;
#endif
	if (max_entropy_shift!=most_5_shift_abs)
	{
		for(unsigned int ic=0; ic<sequences_count; ic++)  //result assign
		//seq-by-seq, including current
		{
			if (!state.motif_present[ic]) continue;
			if (!state.is_complement[ic])
				state[ic]=(original_state[ic]+max_entropy_shift)-most_5_shift_abs;
			else
				state[ic]=(original_state[ic]+most_5_shift_abs)-max_entropy_shift;
		}
		G=Negative_Entropy_with_defined_patterns_per_pattern_position
						(state,max_entropy);
#if __DEBUGLEVEL__ >= 2
		log_stream<<" and G="<<G<<endl;
#endif
		success=have_we_seen_an_allowed_position;
		return 1;
	}
	state=original_state;
	G=Negative_Entropy_with_defined_patterns_per_pattern_position
					(state,max_entropy);
#if __DEBUGLEVEL__ >= 2
	log_stream<<" and G="<<G<<endl;
#endif
	success=have_we_seen_an_allowed_position;
	return 0;
}


unsigned int LookingForAtgcMotifsMultinomialGibbs::
		ShiftAndGapOptimisation
		(
			MarkovChainState & state,
			unsigned int shift_magnitude,
			unsigned short & success,
			//success == 0 if all possibilities are masked
			long double & F,
			long double & best_G
		)
{
	MarkovChainState original_state (state),
									previous_spacer_best_state(sequences_count);

	long double current_spacer_F=0,previous_spacer_F=0,
		previous_spacer_G=0, current_spacer_G=0;
	unsigned int previous_gap=0;
	unsigned int current_gap=0;
	unsigned short is_first_iteration=1;
	unsigned short success_operation;
	success=0;
	do
	{
		previous_gap=current_gap;
		previous_spacer_G=current_spacer_G;
		previous_spacer_F=current_spacer_F;
		previous_spacer_best_state=state;
		state=original_state;
		state.SetSymmetricGap(current_gap);
		//the two previous lines re-initialise state
		if (!is_first_iteration)
		{
			do {state.SetSymmetricGap(++current_gap);}
			while (!state.is_spaced() || //if it is not spaced , repeat
					previous_gap==state.spacer_length());
											//we want the current_gap to change actually
			current_gap=state.spacer_length();
			//it is for safety, virtually it is OK after "while" loop
		}
		if ((state.pattern_length()-current_gap)<shortest_significant_motif)
			break;


		ShiftOptimisation
		(
		 	state,
			shift_magnitude,
			success_operation,
			current_spacer_F,
			current_spacer_G
		);
		success=success||success_operation;
		if (is_first_iteration) {is_first_iteration=0;};
#if __DEBUGLEVEL__ >= 1
		log_stream<<"The gap optimising step. Gap="<<current_gap<<" G="<<current_spacer_G<<" Prev. gap="<<previous_gap<<" Prev. G="<<previous_spacer_G<<endl;
#endif
	} while	(current_gap==0 || previous_spacer_G < current_spacer_G);
	best_G=previous_spacer_G;
	F=previous_spacer_F;
	state=previous_spacer_best_state;
	return (original_state!=state);
}

long double LookingForAtgcMotifsMultinomialGibbs::
	PositionAdjustment
			(
				unsigned long local_steps_made,
				unsigned long adjustments_made,
				unsigned short test_gaps
			)
{
	long double F,G;
	unsigned short success=0;
	if(!be_quiet) log_stream<<"Shift"<<(test_gaps?" and gap ":" ")<<"adjustment (checkpoint after "<<
			local_steps_made<<" steps and "<<adjustments_made<<" regular adjustments made)."<<
			" The state is :"<<endl<<theState;
	//Symbols are OK here.
	G=Negative_Entropy_with_defined_patterns_per_pattern_position(theState,F);
	if(!be_quiet) log_stream<<
			"F is "<<F<<";  information per pattern position (G) is "<<G<<endl
			<<flush;
	//first, we want to find the leftmost and the rightmost possible shift.
	if (!test_gaps)
		last_adjustment_has_changed_state=
				ShiftOptimisation
					(the_state,the_state.pattern_length()-1,success,F,G);
	else
		last_adjustment_has_changed_state=
				ShiftAndGapOptimisation
					(the_state,the_state.pattern_length()-1,success,F,G);

	if (success)
	{
		if (last_adjustment_has_changed_state)
		{
			if(!be_quiet) log_stream<<"The state after adjustment : "<<endl<<theState;
			if(!be_quiet) log_stream<<
					"F is "<<F<<";  information per pattern position is "<<G
					<<flush;
			Symbols->calculate(theState);
		}
		else
		{
			if(!be_quiet) log_stream<<"The adjustment has not change the state."<<endl;
		}
		if(!be_quiet) log_stream<<endl<<endl<<flush;
	}
	else
        throw DumbException("Contradictory state. All shifts are masked.\n");
    return G;
}


long double LookingForAtgcMotifsMultinomialGibbs::
	PositionAndLengthAdjustment
		(
			unsigned int minimal_length,
			unsigned int maximal_length,
			unsigned short test_gaps,
			unsigned long local_steps_made,
			unsigned long adjustments_made
		)
{
	double orig_info;
	double long F;
	unsigned short success=0,success_operation=0;
	if(!be_quiet) log_stream<<"Shift"<<(test_gaps?", gap ":" ")
			<<"and length adjustment (checkpoint after "<<
			local_steps_made<<" steps and "<<adjustments_made<<" regular adjustments made)."<<
			" The state is :"<<endl<<theState;
	orig_info=
			Negative_Entropy_with_defined_patterns_per_pattern_position(theState,F);
	if(!be_quiet) log_stream
			<<"F is "<<F
			<<";  information per pattern position is "<<orig_info<<endl<<flush;
	//first, we want to find the leftmost and the rightmost possible shift.

	if (theState.present_motifs()<2)
	{
		if(!be_quiet) log_stream
				<<"The state is too sparse to do anything."<<endl<<flush;
		last_adjustment_has_changed_state=0;
		return orig_info;
	}
	unsigned int actual_max_length=maximal_length;

	//here, we are looking for maximal possible length
	for (unsigned int i=0; i<sequences_count; i++) //remake!
	//seq-by-seq, including current
	{
		if (!theState.motif_present[i]) continue;
		if (!theState.is_complement[i])
			actual_max_length=
					min(actual_max_length,((int)Sequences[i].size()-theState[i]));
		else
			actual_max_length=
					min(actual_max_length,theState.pattern_length()+theState[i]);
	}
	if (actual_max_length<theState.pattern_length())
	{
        throw DumbException("Contradictory state. Length is more than possible.\n");
    }

    MarkovChainState original_state(theState), optimal_state(theState), state(theState);

    long double max_G=0,curr_G;

	last_adjustment_has_changed_state=0;
	//it's  useful if the minimal_length>maximal_length

	for (
				unsigned int length=minimal_length;
				length<=actual_max_length;
				length++
			)
	{
		//
		//              **********  - start motif
		//              ++++++++++++++ - tested motif
		//
		//              **********
		//                       ++++++++++++++    - most 3'shift
		//
		//
		//              **********
		// ++++++++++++++                          - most 5'shift
		//
		// If we work with complements, the 5' and 3' change their places
		// so,
		int shift_magnitude = min(original_state.pattern_length(),length) - 1;

		state=original_state;

		state.SetLength(length,Sequences,current_caps_mode);

		if (state.present_motifs()<2)break;

		for (unsigned int i=0; i<sequences_count; i++) //remake!
		//we superpose the corresponding motif position after
		//the length change
		{
			if (state.motif_present[i]&&state.is_complement[i])
			{
				state[i]+=theState.pattern_length();
				state[i]-=state.pattern_length();
			}
		}
		//start;

		if (!test_gaps)
			last_adjustment_has_changed_state=
					ShiftOptimisation(state,shift_magnitude,success_operation,F,curr_G);
		else
			last_adjustment_has_changed_state=
					ShiftAndGapOptimisation(state,shift_magnitude,success_operation,F,curr_G);


		if (success_operation) success=1; else continue; // if not the length did not get success, we do not account for it.


		//DEBUG
		//long double G=
		//		Negative_Entropy_with_defined_patterns_per_pattern_position
		//				(state,F);
		//cout<<"best: L="<<length<<"  F="<<F<<" G="<<G<<endl;

		if( length==minimal_length || max_G<curr_G ) //first or better than the best
		{
			optimal_state=state;
			max_G=curr_G;
#if __DEBUGLEVEL__ >= 2
			log_stream<<"+++\nOptimal state assigning. Optimal, the, start (length is "
					<<length<<")\n"<<optimal_state<<state<<original_state<<"+++"<<endl;
#endif
		}
	}//main cycle
	if (success)
	{
		if (optimal_state==original_state)
			last_adjustment_has_changed_state=0;
		else
		{
			last_adjustment_has_changed_state=1;
			the_state=optimal_state;
			Symbols->calculate(the_state);
		}
		if (last_adjustment_has_changed_state)
		{
			if(!be_quiet) log_stream<<"Globally optimal state ( info/pos="<<max_G
					<<" ) is :\n"<<the_state<<endl;
		}
		else
		{
			if(!be_quiet) log_stream<<"Globally optimal state is the same"<<endl;
		}
	}
	else
        throw DumbException("Contradictory state. All shifts are masked.\n");
    /*if (orig_info>max_G) {
		SymbolsCounter * symb1=create_new_counter(Sequences,original_state),
			*symb2=create_new_counter(Sequences,the_state);
		symb1->calculate(original_state);
		symb2->calculate(theState);
		cerr<<orig_info<<">"<<max_G<<endl<<"Orig F="<<symb1->NegativeEntropy()<<
				"  The F="<<symb2->NegativeEntropy()<<endl<<
				original_state<<the_state<<endl<<"BRRRR!!!"<<endl<<flush;
		delete symb1;
		delete symb2;
		exit(1);}*/
    return max_G;
}


void LookingForAtgcMotifsInOneThreadMultinomialGibbs::
	calculate_position_weights
		(
			const MarkovChainState & state,
			unsigned int sequence_no,
			unsigned int motif_positions,
			SymbolsCounter & symbols,
			//use this SymbolsRate for count position weights,
			vector<double> & y
			//place the results here.
			//the vector it to be at least
			//last_possible_position-first_possible_position+1
		)
		//symbols is already counted by	calculate() and exclude_sequence(...)
		//
		//the result we wait for is:
		//y[0] is posterior for motif absence in sequence_no
		//y[1..motif_positions] are posteriors for
		//motif presence in the position pos-1 (0-based)
{
#if __DEBUGLEVEL__ >= 1
	SymbolsCounter * ac=create_new_counter(Sequences,state);
	ac->calculate(state);
	ac->exclude_sequence(state,sequence_no);

	if (!(*ac == symbols)) log_stream<<"Given and paramter-based counters discrepancy in calc_weights 1 thread)\n"<<
		"\ncurrent sequence is "<<sequence_no<<
				"\nstate is:"<<state<<
				"\ngiven:\n"<<symbols<<"\ncalculated:\n"<<*ac<<"\n"<<flush;
	delete ac;
#endif
	vector <double> a((letters_in_alphabet+1)*state.pattern_length());
	//it is matrix of q(pos,symbol)/p(symbol)
	//maybe later it will be split to a and p,
	//but now it is  not necessary
	//a(pos,symbol) is a[(symbol-1)*pattern_length+pos]
	//pos in pattern is 0-based, symbol is 1-based.
	for (unsigned short j=1;j<=letters_in_alphabet;j++)
		for (unsigned int i=0;i<state.pattern_length();	i++)
			a[(j)*state.pattern_length()+i]=
					symbols.foreground_probability(i,j)/
					symbols.background_probability(j); //i is pos, j is symbol

	for (unsigned int i=0;i<state.pattern_length();	i++)
		a[i]=0.;

	for (unsigned int p=1;p<=motif_positions;p++)
	//pattern positions, one-by-one
	{
		//
		double A=1.;
		if (current_caps_mode==off)
		{
			for (unsigned int i=0;i<state.pattern_length();i++)
			{
				A*=a[( Sequences[sequence_no][(p-1)+i] )*state.pattern_length()+i];
				//q_p_j/p_j;
				//if we meet 0 (mask) a[0+i] is 0 for any i in the cycle band and
				//so we get 0. Mask is coded as 0 in the sequence.
			}
			y[p]=A;
			//it is posterior(motif_in_this_position)/posterior(no_motif)
		}
		else
		{
			unsigned int caps_counter=0;
			for (unsigned int i=0;i<state.pattern_length();i++)
				caps_counter+=Sequences.caps[sequence_no][(p-1)+i];
			if (
						( current_caps_mode==one && caps_counter)
						||
						(current_caps_mode==all && caps_counter==state.pattern_length())
					)
			{
				for (unsigned int i=0;i<state.pattern_length();i++)
					A*=a[( Sequences[sequence_no][(p-1)+i] )*state.pattern_length()+i];
					//q_p_j/p_j;
					//if we meet 0 (mask) a[0+i] is 0 for any i in the cycle band and
					//so we get 0. Mask is coded as 0 in the sequence.
				y[p]=A;
				//it is posterior(motif_in_this_position)/posterior(no_motif)
				//if there were no caps, the probability is 0
			}
			else
				y[p]=0.;
		}
	}
	y[0] = motif_positions*motif_absence_prior/(1-motif_absence_prior);
	//it is prior(no_motif)/prior(motif_in_definite_position)
	//that's it
}

void LookingForAtgcMotifsInTwoThreadsMultinomialGibbs::
	calculate_position_weights
		(
			const MarkovChainState & state,
			unsigned int sequence_no,
			unsigned int motif_positions,
			SymbolsCounter & symbols,
			//use this SymbolsRate for count weights,
			vector<double> & y
			//place the results here.
			//the vector it to be at least
			//last_possible_position-first_possible_position+1
		)
		//symbols is already counted by	calculate() and exclude_sequence(...)
		//
		//the result we wait for is:
		//y[0] is posterior for motif absence in sequence_no read in direct way
		//y[1..motif_positions] are posteriors for
		//motif presence in the position pos-1 (0-based) in the same way
		//y[motif_positions+1..2*motif_positions are posteriors of
		//motif present in position pos-motif_positions-1 (0-based) read as
		//complement
		//y[2*motif_positions+1] is posterior of complementary sequence
		//and motif absence.
{
// the things remain the same between motif positions and motif absence.
// all the positions which suppose to read the sequence as complement
// are to be multiplied by       mult (P(compl(a))/p(a)
//                        0..length-1

#if __DEBUGLEVEL__ >= 1
	SymbolsCounter * ac=create_new_counter(Sequences,state);
	ac->calculate(state);
	ac->exclude_sequence(state,sequence_no);
	if (!(*ac == symbols)) log_stream<<"Given and paramter-based counters discrepancy in calc_weights 2 thread\n"<<
		"current sequence is "<<sequence_no<<
				"\nstate is:"<<state<<
				"\ngiven:\n"<<symbols<<"\ncalculated:\n"<<*ac<<"\n"<<flush;
	delete ac;
#endif
	//cout<<"%"<<current_caps_mode<<endl<<flush;
	vector <double> a((letters_in_alphabet+1)*state.pattern_length());
	vector <double> p(letters_in_alphabet+1);
	//it is matrix of q(pos,symbol)/p(symbol)
	//maybe later it will be split to a and p,
	//but now it is  not necessary
	//a(pos,symbol) is a[(symbol-1)*pattern_length+pos]
	//pos in pattern is 0-based, symbol is 1-based.
	for (unsigned int i=0;i<state.pattern_length();	i++)
		a[i]=0.;
	p[0]=1.;
	for (unsigned short j=1;j<=letters_in_alphabet;j++)
	{
		p[j]=symbols.background_probability(j);
		for (unsigned int i=0;i<state.pattern_length();	i++)
			a[j*state.pattern_length()+i]=
					symbols.foreground_probability(i,j)/p[j];
	}
	p[1]=p[2]/p[1];
	p[2]=1./p[1];
	p[3]=p[4]/p[3];
	p[4]=1./p[3];
	//now, every p is p(compl(a))/p(a)
	//we need it only for mult(all p(compl(b)/p(b))
	//here, we count the coefficient we are to multiply all
	//posteriors dealing with complement.
	double complement_multiplier=1.;
	for (unsigned int i=0;i<Sequences[sequence_no].size();i++)
		complement_multiplier*=p[Sequences[sequence_no][i]];
	for (unsigned int p=1;p<=motif_positions;p++)
	//pattern positions, one-by-one
	{
		//
		double A=1.,AC=1.;
		if (current_caps_mode==off)
		{
			for (unsigned int i=0;i<state.pattern_length();i++)
			{
				A*=a[( Sequences[sequence_no][(p-1)+i] )*state.pattern_length()+i];
				//q_p_j/p_j;
				//if we meet 0 (mask) a[0+i] is 0 for any i in the cycle band and
				//so we get 0. Mask is coded as 0 in the sequence.
				AC*=a[
						(
							Atgc::complement
								(
									Sequences[sequence_no][(p-1)+i]
								)
						)*state.pattern_length()+
						state.pattern_length()-i-1
					];
			}
			y[p]=A;
			y[p+motif_positions]=AC*complement_multiplier;
			//it is posterior(motif_in_this_position)/posterior(no_motif)
		}
		else
		{
			unsigned int caps_counter=0;
			for (unsigned int i=0;i<state.pattern_length();i++)
				caps_counter+=Sequences.caps[sequence_no][(p-1)+i];
			if (
						( current_caps_mode==one && caps_counter)
						||
						(current_caps_mode==all && caps_counter==state.pattern_length())
					)
			{
				for (unsigned int i=0;i<state.pattern_length();i++)
				{
					A*=a[( Sequences[sequence_no][(p-1)+i] )*state.pattern_length()+i];
					//q_p_j/p_j;
					//if we meet 0 (mask) a[0+i] is 0 for any i in the cycle band and
					//so we get 0. Mask is coded as 0 in the sequence.
					AC*=a[
							(
								Atgc::complement
									(
										Sequences[sequence_no][(p-1)+i]
									)
							)*state.pattern_length()+
							state.pattern_length()-i-1
						];
				}
				y[p]=A;
				//it is posterior(motif_in_this_position)/posterior(no_motif)
				//if there were no caps, the probability is 0
				y[p+motif_positions]=AC*complement_multiplier;
			}
			else
			{
				y[p]=0.;
				y[p+motif_positions]=0.;
			}
		}
	}
	y[0] = motif_positions*motif_absence_prior/(1-motif_absence_prior);
	//it is prior(no_motif)/prior(motif_in_definite_position)
	y[2*motif_positions+1] =
			complement_multiplier*motif_positions*motif_absence_prior/
			(1-motif_absence_prior);
	//it is poterior for absence - and - comlement chain
	//that's it
}

void LookingForAtgcPalindromesMultinomialGibbs::
	calculate_position_weights
		(
			const MarkovChainState & state,
			unsigned int sequence_no,
			unsigned int motif_positions,
			SymbolsCounter & symbols,
			//use this SymbolsCounter for count the weights,
			vector<double> & y
			//place the results here.
			//the vector it to be at least
			//last_possible_position-first_possible_position+1
		)
		//symbols is already counted by	calculate() and exclude_sequence(...)
		//
		//the result we wait for is:
		//y[0] is posterior for motif absence in sequence_no
		//y[1..motif_positions] are posteriors for
		//motif presence in the position pos-1 (0-based)
{

#if __DEBUGLEVEL__ >= 1
	SymbolsCounter *ac=create_new_counter(Sequences,state);
	ac->calculate(state);
	ac->exclude_sequence(state,sequence_no);

	if (!(*ac == symbols)) log_stream<<"Given and paramter-based counters discrepancy in calc_weights palindromes!!!\n"<<
		"\ncurrent sequence is "<<sequence_no<<
				"\nstate is:"<<state<<
				"\ngiven:\n"<<symbols<<"\ncalculated:\n"<<*ac<<"\n"<<flush;
#endif
	vector <double> a((letters_in_alphabet+1)*state.pattern_length());
	//it is matrix of q(pos,symbol)/p(symbol)
	//maybe later it will be split to a and p,
	//but now it is  not necessary
	//a(pos,symbol) is a[(symbol-1)*pattern_length+pos]
	//pos in pattern is 0-based, symbol is 1-based.
	for (unsigned short j=1;j<=letters_in_alphabet;j++)
		for (unsigned int i=0;i<state.pattern_length();	i++)
				a[j*state.pattern_length()+i]=
					symbols.foreground_probability(i,j)/
					symbols.background_probability(j); //i is pos, j is symbol

	for (unsigned int i=0;i<state.pattern_length();	i++)
		a[i]=0.;

	for (unsigned int p=1;p<=motif_positions;p++)
	//pattern positions, one-by-one
	{
		//
		double A=1.;
		if (current_caps_mode==off)
		{
			for (unsigned int i=0;i<state.pattern_length();i++)
			{
				A*=a[( Sequences[sequence_no][(p-1)+i] )*state.pattern_length()+i];
				//q_p_j/p_j;
				//if we meet 0 (mask) a[0+i] is 0 for any i in the cycle band and
				//so we get 0. Mask is coded as 0 in the sequence.
			}
			y[p]=A;
			//it is posterior(motif_in_this_position)/posterior(no_motif)
		}
		else
		{
			unsigned int caps_counter=0;
			for (unsigned int i=0;i<state.pattern_length();i++)
				caps_counter+=Sequences.caps[sequence_no][(p-1)+i];
			if (
						( current_caps_mode==one && caps_counter)
						||
						(current_caps_mode==all && caps_counter==state.pattern_length())
					)
			{
				for (unsigned int i=0;i<state.pattern_length();i++)
					A*=a[( Sequences[sequence_no][(p-1)+i] )*state.pattern_length()+i];
					//q_p_j/p_j;
					//if we meet 0 (mask) a[0+i] is 0 for any i in the cycle band and
					//so we get 0. Mask is coded as 0 in the sequence.
				y[p]=A;
				//it is posterior(motif_in_this_position)/posterior(no_motif)
				//if there were no caps, the probability is 0
			}
			else
				y[p]=0.;
			//it is posterior(motif_in_this_position)/posterior(no_motif)
		}
	}
	y[0] = motif_positions*motif_absence_prior/(1-motif_absence_prior);
	//it is prior(no_motif)/prior(motif_in_definite_position)
	//that's it
}


unsigned int LookingForAtgcMotifsMultinomialGibbs::
	get_int_draw_from_a_ditribution
//it generates outcome from a distribution represented
//by first_possible_position..last_possible_position
//in Y
		(
			unsigned int least_possible_draw,
			unsigned int most_possible_draw,
			vector<double> y
			//get the distribution from here.
			//the vector it to be at least
			//last_possible_position-first_possible_position+1
		)
{
	double draw=uni();
	double sum=0.;
	//scaling the draw
	long double weights_sum=0.;
	for (unsigned int ij=least_possible_draw;ij<=most_possible_draw;ij++)
		weights_sum+=y[ij];
	draw*=weights_sum;
	//scaled
	unsigned int i=least_possible_draw;
	while (sum + y[i] < draw && i<=most_possible_draw) sum+=y[i++];
	//so, the probability to dtop on position i is proportinal to
	//y[i].
  //
	//|                  |     |          |     |
	//      y[0]           y[1]   y[2]       y[4]
	//                                  +
	//                                 draw
	//
	//y[0]+y[1] < draw ; y[0]+y[1]+y[2] >draw => i=2
	//
	if (i==most_possible_draw+1)
	{
        throw DumbException("Trying to make draw more than the most possible.\n");
    }
    return i;
}

double long LookingForAtgcMotifsMultinomialGibbs::
	Negative_Entropy_with_defined_patterns_per_pattern_position
		(
			const MarkovChainState & state,
			double long & NegativeEntropy
		)
{
	//When counting the space distribution information gain,
	//we use a simplification:
	//
	//(sum( (y_i/w) log(y_i/w) ), where w=sum (y_i) equals
	//    (1/w)( sum (y_i (log(y_i)-log(w) ) =
	//=   (1/w) sum (y_i log(y_i)) - (1/w log (w) sum (y_i)) =
	//=   (1/w) sum (y_i log(y_i)) - log(w)
	//
	double long result;
	SymbolsCounter * sym_w = create_new_counter(Sequences,state);
	sym_w->calculate(state);
	double long Position_Information=0.;
	double long position_information_for_the_sequence;
	double long prob_coeff;
	double long weights_sum;
//	unsigned int informative_part_length=
//			state.pattern_length()-(state.is_spaced()?
//					(state.spacer_3_end()-state.spacer_5_end()+1):0);

	for (unsigned i=0;i<sequences_count;i++) //sequence - by - sequence
	{
		if (state.motif_was_too_long[i]) continue;
		position_information_for_the_sequence=0.;
		unsigned int motif_positions=Sequences[i].size()-
				state.pattern_length()+1;
		unsigned int possible_motif_positions;
		sym_w->exclude_sequence(state,i);
		unsigned int leftmost_masked_position;
		//taking masked positions out of account
		//
		if (mode==two_threads)
		{
			possible_motif_positions=2*motif_positions;

			//calculating masked positions for direct strand
			leftmost_masked_position=Sequences[i].size()+1;
			for(unsigned int p=0;p<state.pattern_length();p++)
				if(Sequences.mask[i][p]&1u) leftmost_masked_position=p;

			for(unsigned int j=0;j<motif_positions;j++)
			{
				if(Sequences.mask[i][j+state.pattern_length()-1]&1u)
					leftmost_masked_position=j+state.pattern_length()-1;
				if (leftmost_masked_position>=j &&
						leftmost_masked_position<=j+state.pattern_length()-1)
					possible_motif_positions--;
			}
			//calculating masked positions for opposite strand
			//(we go in usual direction to simplify the life)
			leftmost_masked_position=Sequences[i].size()+1;
			for(unsigned int p=0;p<state.pattern_length();p++)
				if(Sequences.mask[i][p]&2u) leftmost_masked_position=p;

			for(unsigned int j=0;j<motif_positions;j++)
			{
				if(Sequences.mask[i][j+state.pattern_length()-1]&2u)
					leftmost_masked_position=j+state.pattern_length()-1;
				if (leftmost_masked_position>=j &&
						leftmost_masked_position<=j+state.pattern_length()-1)
					possible_motif_positions--;
			}

			calculate_position_weights
					(
						state,
						i,
						motif_positions,
						*sym_w, //we use private symbols counter not to destroy Symbols
						Y
					);
			//count weights_sum
			weights_sum=0;
			for(unsigned int j=0;j<=2*motif_positions+1;j++)
				weights_sum+=Y[j];
			//
			for(unsigned int j=1;j<=2*motif_positions;j++)
			{
				prob_coeff=Y[j]*possible_motif_positions/(1-motif_absence_prior);
				if (prob_coeff>EPSILON)
					position_information_for_the_sequence-=Y[j]*log_2(prob_coeff);
			}
			if (motif_absence_prior>EPSILON)
			{
				prob_coeff=Y[0]*2./motif_absence_prior;
				if (prob_coeff>EPSILON)
					position_information_for_the_sequence-=Y[0]*log_2(prob_coeff);
				prob_coeff=Y[2*motif_positions+1]*2./motif_absence_prior;
				if (prob_coeff>EPSILON)
					position_information_for_the_sequence-=
												Y[2*motif_positions+1]*log_2(prob_coeff);
			}
			position_information_for_the_sequence/=weights_sum;
			position_information_for_the_sequence+=log_2(weights_sum);
			//
			//because -SUM_i(y_i/sum_j(y_j)log(y_i/(p_i*sum_j(y_j))) =
			// =log(sum_j(y_j)-(1/sum_j(y_j))*SUM_i(y_i*log(y_i/p_i))
			//
		}
		else //not two_threads
		{
			possible_motif_positions=motif_positions;
			//calculating masked positions for direct strand
			leftmost_masked_position=Sequences[i].size()+1;
			for(unsigned int p=0;p<state.pattern_length();p++)
				if(Sequences.mask[i][p]&1u) leftmost_masked_position=p;

			for(unsigned int j=0;j<motif_positions;j++)
			{
				if(Sequences.mask[i][j+state.pattern_length()-1]&1u)
					leftmost_masked_position=j+state.pattern_length()-1;
				if (leftmost_masked_position>=j &&
						leftmost_masked_position<=j+state.pattern_length()-1)
					possible_motif_positions--;
			}
			calculate_position_weights
					(
						state,
						i,
						motif_positions,
						*sym_w, //we use private symbols counter not to destroy Symbols
						Y
					);
			//count weights_sum
			weights_sum=0;
			for(unsigned int j=0;j<=motif_positions;j++)
				weights_sum+=Y[j];
			//
			for(unsigned int j=1;j<=motif_positions;j++)
			{
				prob_coeff=Y[j]*possible_motif_positions/(1-motif_absence_prior);
				if (prob_coeff>EPSILON)
					position_information_for_the_sequence-=Y[j]*log_2(prob_coeff);
			}
			if (motif_absence_prior>EPSILON)
			{
				prob_coeff=Y[0]/motif_absence_prior;
				if (prob_coeff>EPSILON)
					position_information_for_the_sequence-=Y[0]*log_2(prob_coeff);
			}
			position_information_for_the_sequence/=weights_sum;
			position_information_for_the_sequence+=log_2(weights_sum);
			//
			//because -SUM_i(y_i/sum_j(y_j)log(y_i/(p_i*sum_j(y_j))) =
			// =log(sum_j(y_j)-(1/sum_j(y_j))*SUM_i(y_i*log(y_i/p_i))
			//
		}
		sym_w->include_sequence(state,i);
		Position_Information+=position_information_for_the_sequence;
	}

  NegativeEntropy=sym_w->NegativeEntropy();
//	cout<<"Position information: "<<Position_Information<<
//			"  Structural information: "<<NegativeEntropy<<endl<<flush;
	result=(NegativeEntropy+Position_Information)/
			((letters_in_alphabet-1)*state.informative_part_length()
			 *sym_w->motifs());
	delete sym_w;
#if __DEBUGLEVEL__>=2
	log_stream<<"\n$$$$\n  The info per freedom degree is "<<result<<
			"=("<<NegativeEntropy<<"+"<<Position_Information<<")/"<<
						((letters_in_alphabet-1)*sym_w->motifs()*
						 			 state.informative_part_length())<<" for state\n    "<<state<<"\n$$$$\n"<<endl;
#endif
	return result;
}


void LookingForAtgcMotifsInTwoThreadsMultinomialGibbs::
	BalanceComplements()
//nothing mystical:
//we want to awoid the sutuation when we read more then half
//sequences as complements :
//it is more polite to complement all the sequences is the case
{
	if(!be_quiet) log_stream<<"Balancing complements..."<<endl<<flush;
	if (
			(unsigned int)count
			(
				theState.is_complement.begin(),
				theState.is_complement.end(),
				1
			)*2 >= sequences_count
		) //more than half of seqs are complements
	{
		for (unsigned int seq=0;seq<sequences_count;seq++)
		{
			the_state.is_complement[seq]=!the_state.is_complement[seq];
			(static_cast<AtgcSymbolsCounter*>(Symbols))->
					current_sequence_was_complemented(seq);
		}
		Symbols->calculate(theState);
	}
}

void LookingForAtgcMotifsInTwoThreadsMultinomialGibbs::NormaliseState()
{
	BalanceComplements();
}
void LookingForAtgcPalindromesMultinomialGibbs::NormaliseState(){};
void LookingForAtgcMotifsInOneThreadMultinomialGibbs::NormaliseState(){};

