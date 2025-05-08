/****************************************************************************\
SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2024
$Id$
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
#include <algorithm>
#include <functional>
//#include <filesystem>
using namespace std;

#include "Sequences.hpp"
#include "MarkovChainState.hpp"
#include "ResultsSet.hpp"
#include "MCMC.hpp"
#include "Logger.hpp"
#include "Diagnostics.hpp"

extern char * core_git_id;

const char * version_name="SeSiMCMC 4.36.1, fur Alyza";

extern "C"
{
	#include "Random.h"
	#include "confread.h"
}

unsigned int fake_data_positions[20]={84,81,78,75,72,69,66,63,60,50,
													47,44,41,38,35,32,29,26,23,20};

void createFakeExperimentalBunch(SequencesPile &sp,unsigned int be_quiet=0,
																	ostream * log_file_ptr=&cerr,
																	double noise_prob=0);

void createFakeExperimentalGappedBunch(SequencesPile &sp,unsigned int be_quiet=0,
																	ostream * log_file_ptr=&cerr,
																	double noise_prob=0);

void help()
{
	cout<<
	"SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2021."<<endl<<
	endl<<
	"Usage : SeSiMCMC [options] [(--config-file|-f) config-file] < (FastA input)."<<endl<<
	endl<<
	"SeSiMCMC < (FastA input) or"<<endl<<""
	"SeSiMCMC -i FastA_input_file is a good how-do-you-do."<<endl<<
	endl<<
	"To get a complete configuration file, say -w in addition to other options:"<<endl<<
	"the configuration file of the task will be outputed instead of results."<<endl<<
	endl<<
	"Some options that affect only output form, like --html or --cgi can be provided"<<endl<<
	"only via command-line keys."<<endl<<
	endl<<
	"For options and configuration tags information, refer to readme file."<<endl<<
	endl<<
	endl;
}

void version()
{
	char core_id_string[120]; //too much but who knows
	strncpy(core_id_string, core_git_id, 119);
	core_id_string[119]='\0';
	char * core_version;
	strtok(core_id_string," ");
	core_version=strtok(NULL," ");

	cout<<version_name<<endl<<
	"core (MCMC.cpp) git id is "<<core_version<<endl;
}

struct config_parameters
{
	config_parameters():background()
	{
		mode=two_threads;
		smode=no;
		ctmode=information;
		motif_absence_prior=1.1;
		adjust_motif_length=1;
		spaced_motif=0;
		initial_caps=off;
		annealing_caps=off;
		refinement_caps=off;
		adaptive_pseudocouts=0;
		common_background=0;
		motif_length=0;
		maximal_motif_length=0;
		minimal_motif_length=0;
		retrieve_mode=undef_retrieve;
		retrieve_threshold=0x10000;//default will be 1 later
		trim_edges=0;
		dumb_trim_edges=0;
		local_step_cycles_between_adjustments=0;
		adjustments_during_annealing=1;
		cycles_with_minor_change_to_say_cold=0; //insensible
		minor_change_level=-1.; //this, too
		annealings_number_maximum_is_to_be_global_for=10;
		steps_number_maximum_is_to_be_global_for=0;
		slow_optimisation=0;
		chain_fails_after=3;
		chains_to_try=3;
		cycles_per_annealing_attempt=0;
		annealing_attempts=5;
		seed1=1248312;
		seed2=7436111;
/*****************************************************************************/
/*here, we put the defaults for parameters which do not depend on data          */
/*and that are given as parameters to routines that find best motif             */
/*****************************************************************************/
//	annealings_number_maximum_is_to_be_global_for=10;
//	cycles_with_minor_change_to_say_cold=3;
//	minor_change_level=0.5;
//	int cycles_with_minor_change_to_say_cold=5;
//	minor_change_level=0.7;
//	cycles_with_minor_change_to_say_cold=0; //insensible
//	minor_change_level=-1.; //this, too
//	ctmode=informational;
//	chain_fails_after=3;
//  chains_to_try=3;
//	seed1=1248312,seed2=7436111;
/*****************************************************************************/
/*here, we put the the fake defaults which depends on sequence data and their      */
/*default values are calculated in sampler constructor                                             */
/*****************************************************************************/
//	motif_length=0;
//	maximal_motif_length=0;
//	minimal_motif_length=0;
//	local_step_cycles_between_adjustments=0;
//	motif_absence_prior=1.1;
//	cycles_per_annealing_attempt=0;
/*****************************************************************************/
	}
	complement_mode mode;
	symmetry_mode smode;
	change_test_mode ctmode;
	double motif_absence_prior;
	unsigned int adjust_motif_length;
	unsigned int spaced_motif;
	caps_mode initial_caps;
	caps_mode annealing_caps;
	caps_mode refinement_caps;
	unsigned int adaptive_pseudocouts;
	unsigned int common_background;
	unsigned int motif_length;
	unsigned int maximal_motif_length;
	unsigned int minimal_motif_length;
	retrieve_mode_type retrieve_mode;//default smart
	unsigned int retrieve_threshold;//default 0 or 1
	unsigned int trim_edges; //default 0
	unsigned int dumb_trim_edges; //default 0
	unsigned long local_step_cycles_between_adjustments;
	unsigned int adjustments_during_annealing;
	unsigned int cycles_with_minor_change_to_say_cold; //insensible
	double minor_change_level; //this, too
	unsigned int annealings_number_maximum_is_to_be_global_for;
	unsigned long steps_number_maximum_is_to_be_global_for;
	unsigned int slow_optimisation;
	unsigned int chain_fails_after;
	unsigned int chains_to_try;
	unsigned long cycles_per_annealing_attempt;
	unsigned int annealing_attempts;
	unsigned long seed1,seed2;
	vector<double> background;

	string & cipher_caps_mode()
	{
		string & mode = * new string;
		if (initial_caps==annealing_caps and annealing_caps==refinement_caps)
		{
			switch (initial_caps)
			{
			case off: mode="off";return mode;
			case one: mode="one";return mode;
			case all: mode="all";return mode;
			default: mode="unknown_caps_mode";return mode;
			}
		}
		//if we are here, they differ
		switch (initial_caps)
		{
		case off: mode="off-";break;
		case one: mode="one-";break;
		case all: mode="all-";break;
		default: mode="unknown_caps_mode";return mode;
		}
		switch (annealing_caps)
		{
		case off: mode+="off-";break;
		case one: mode+="one-";break;
		case all: mode+="all-";break;
		default: mode="unknown_caps_mode";return mode;
		}
		switch (refinement_caps)
		{
		case off: mode+="off";break;
		case one: mode+="one";break;
		case all: mode+="all";break;
		default: mode="unknown_caps_mode";return mode;
		}
		return mode;
	}

	unsigned int decipher_caps_mode(const string mode)
	//0 is OK
	//1 is format error
	//2 is wrong mode (e.g with increasing restrictivity)
	{
		caps_mode i_c,a_c,r_c;
		if (mode.length()==3)
		{
			if (mode=="off")
				i_c=r_c=a_c=off;
			else if (mode=="one")
				i_c=r_c=a_c=one;
			else if (mode=="all")
				i_c=r_c=a_c=all;
			else return 1;
		}
		else if (mode.length()==11)
		{
			string mode_i=mode.substr(0,3);
			string mode_a=mode.substr(4,3);
			string mode_r=mode.substr(8,3);
			if (mode_i=="off")
				i_c=off;
			else if (mode_i=="one")
				i_c=one;
			else if (mode_i=="all")
				i_c=all;
			else return 1;
			if (mode_a=="off")
				a_c=off;
			else if (mode_a=="one")
				a_c=one;
			else if (mode_a=="all")
				a_c=all;
			else return 1;
			if (mode_r=="off")
				r_c=off;
			else if (mode_r=="one")
				r_c=one;
			else if (mode_r=="all")
				r_c=all;
			else return 1;
			if (r_c>a_c || a_c>i_c) return 2;
		}
		else return 1;
		initial_caps=i_c;
		annealing_caps=a_c;
		refinement_caps=r_c;
		return 0;
	}


	ostream & print(ostream & o, string prefix_string="", string posfix_string="",
			                         string each_string_prefix="")
			//prefix_string and posfix_string are to contain "\n" if they are strings
	{
		string & caps_mode(cipher_caps_mode());
		o<<prefix_string;
		o<<each_string_prefix<<"mode=";
		switch (mode)
		{
		case one_thread: o<<"one_strand"<<endl;break;
		case two_threads: o<<"two_strands"<<endl;break;
		case unknown_mode:o<<"unknown"<<endl;
		}
		o<<each_string_prefix<<"symmetry=";
		switch (smode)
		{
		case no: o<<"no"<<endl;break;
		case repeats: o<<"repeats"<<endl;break;
		case palindromes: o<<"palindromes"<<endl;break;
		case unknown_smode:o<<"unknown"<<endl;
		}
		o<<each_string_prefix<<"motif_absence_prior="<<motif_absence_prior<<endl;
		o<<each_string_prefix<<"adjust_motif_length="<<(adjust_motif_length?"yes":"no")<<endl;
		o<<each_string_prefix<<"spaced_motif="<<(spaced_motif?"yes":"no")<<endl;

		o<<each_string_prefix<<"common_background="<<(common_background?"yes":"no")<<endl;
		o<<each_string_prefix<<"adaptive_pseudocouts="<<(adaptive_pseudocouts?"yes":"no")<<endl;
		o<<each_string_prefix<<"motif_length="<<motif_length<<endl;
		o<<each_string_prefix<<"maximal_motif_length="<<maximal_motif_length<<endl;
		o<<each_string_prefix<<"minimal_motif_length="<<minimal_motif_length<<endl;
		switch (retrieve_mode)
		{
			case no_retrieve:
				o<<each_string_prefix<<"retrieve_other_sites=no"<<endl;
				break;
			case yes_retrieve:
				o<<each_string_prefix<<"retrieve_other_sites=yes"<<endl;
				if (retrieve_threshold != 0x10000)
					o<<each_string_prefix<<"retrieve_other_sites_threshold="<<retrieve_threshold<<endl;
				break;
			case smart_retrieve:
				o<<each_string_prefix<<"retrieve_other_sites=smart"<<endl;
				if (retrieve_threshold != 0x10000)
					o<<each_string_prefix<<"retrieve_other_sites_threshold="<<retrieve_threshold<<endl;
				break;
			case pure_smart_retrieve:
				o<<each_string_prefix<<"retrieve_other_sites=pure_smart"<<endl;
				break;
			default:;
		}
		o<<each_string_prefix<<"trim_edges="<<(trim_edges?"yes":"no")<<endl;
		o<<each_string_prefix<<"dumb_trim_edges="<<(dumb_trim_edges?"yes":"no")<<endl;
		o<<each_string_prefix<<"caps="<<caps_mode<<endl;
		o<<each_string_prefix<<"local_step_cycles_between_adjustments="<<
				local_step_cycles_between_adjustments<<endl;
		o<<each_string_prefix<<"adjustments_during_annealing="<<
				(adjustments_during_annealing?"yes":"no")<<endl;
		o<<each_string_prefix<<"annealing_criterion=";
		switch (ctmode)
		{
		case site_positions: o<<"site_positions";break;
		case information: o<<"information";break;
		case correlation: o<<"correlation";break;
		case unknown_ctmode:;
		}
		o<<endl;
		o<<each_string_prefix<<"cycles_with_minor_change_to_say_cold="<<
				cycles_with_minor_change_to_say_cold<<endl;
		o<<each_string_prefix<<"minor_change_level="<<minor_change_level<<endl;
		o<<each_string_prefix<<"annealings_number_maximum_is_to_be_global_for="<<
				annealings_number_maximum_is_to_be_global_for<<endl;
		o<<each_string_prefix<<"steps_number_maximum_is_to_be_global_for="<<
				 steps_number_maximum_is_to_be_global_for<<endl;
		o<<each_string_prefix<<"slow_optimisation="<<(slow_optimisation?"yes":"no")<<endl;
		o<<each_string_prefix<<"chain_fails_after="<<chain_fails_after<<endl;
		o<<each_string_prefix<<"chains_to_try="<<chains_to_try<<endl;
		o<<each_string_prefix<<"cycles_per_annealing_attempt="<<cycles_per_annealing_attempt<<endl;
		o<<each_string_prefix<<"annealing_attempts="<<annealing_attempts<<endl;
		o<<each_string_prefix<<"random_seed_1="<<seed1<<endl;
		o<<each_string_prefix<<"random_seed_2="<<seed2<<endl;
		if (background.size()==4)
		{
			o<<each_string_prefix<<"background_A="<<background[0]<<endl;
			o<<each_string_prefix<<"background_T="<<background[1]<<endl;
			o<<each_string_prefix<<"background_G="<<background[2]<<endl;
			o<<each_string_prefix<<"background_C="<<background[3]<<endl;
		}
		o<<posfix_string;
		return o;
	}
};

/*****************************************************************************
 * Here we put dictionary of old -> new switches/config tags because we do
 * not want to change the varaibles inside the code.
 * --one-thread  -> --one-strand
 * --two-threads -> --two-strands
 * --no-gaps     -> --no-spacers
 * --not-gapped  -> --not-spaced
 * --show-gapped-bases
 *               -> --show-spaced-bases
 * --sgb         -> --ssb
 * --hide-gapped-bases
 *               -> --hide-spaced-bases
 * --not-show-gapped-bases
 *               -> --not-show-spaced-bases
 * --hgb         -> --hsb
 *****************************************************************************/

int main(int argc, char ** argv)
{
    try {
        config_parameters config;
        char mode_string[100];
        LookingForAtgcMotifsMultinomialGibbs *sampler = NULL;

        /*****************************************************************************/
        /*here, we put the defaults for all things, except those which depend on data*/
        /*****************************************************************************/
        unsigned int do_nothing = 0;
        unsigned long time_limit = 1000;
        unsigned int time_limit_is_default = 1;
        unsigned int be_quiet = 1;
        unsigned int be_not_quiet = 0;
        //	unsigned long steps_number_maximum_is_to_be_global_for=0;
        //	unsigned int adjust_motif_length=1;
        //	unsigned int spaced_motif=0;
        //	complement_mode mode=two_threads;
        //	symmetry_mode smode=no;
        unsigned short generate_fake_data = 0;
        //	unsigned int slow_optimisation=0;
        //	unsigned int common_background=0;
        unsigned int given_background_counters = 0;
        //	vector<double> background;
        double background_A;
        double background_T;
        double background_G;
        double background_C;
        unsigned int flanks_length = 5;         //default flanks
        unsigned short show_letters_in_gap = 1; //default flanks
        unsigned int background_A_set = 0;
        unsigned int background_T_set = 0;
        unsigned int background_G_set = 0;
        unsigned int background_C_set = 0;
        unsigned int output_xml = 0;
        unsigned int output_html = 0;
        unsigned int output_cgi = 0;
        unsigned int output_weight_tables = 0;
        unsigned int not_output_weight_tables = 0;
        //ibis format pwm output
        unsigned int ibis_output = 0;
        string TF_name = "noname";
        string motif_suffix = "noname";
        //ibis format optins end
        //	unsigned int adjustments_during_annealing=1;
        string InputFileName = "";
        string FakeInputFileName = "";
        string OutputFileName = "";
        string OutputFastAFileName = "";
        double masked_part_on_output = 0.5;
        unsigned short OutputFileOpenedOK = 0;
        unsigned short OutputFastAFileOpenedOK = 0;
        string logfile = "/dev/null";
        string ip = "unknown";
        string id = "undefined";
        istream *fastaptr = &cin;
        unsigned int errlevel = 0;
        Diagnostics diagnostics;

        unsigned int output_unreliable_result = 0;

        AcgtLetterOrder order;

        ostream *log_file_ptr = &cerr;
        ostream *out_ptr = &cout;
        ofstream *outf_ptr = NULL;
        //if we output to a file, we make all the file-related operations
        //using outf_ptr and we ouutput to out_ptr;
        //actually, they refer to the same object.
        //It is neccessary to use the same output(something,out_ptr)
        //interface both for file-based and stream-based output
        ofstream *outfasta_ptr = NULL;

        //the set is necessary for slow mode results output
        set<pair<unsigned int, double> > G_values;

        const char *cfilenameptr = NULL;
        FILE *cfile = NULL;
        unsigned short config_file_given = 0;
        unsigned short input_file_name_rewritten = 0;
        unsigned int arguments_reflected = 0;
        //parsing arguments
        vector<string> args;
        vector<string>::iterator op;
        const char *value;
        char *ptr;

        unsigned int forbid_permutations = 0; //test workaround

        for (int argn = 1; argn < argc; argn++)
            args.push_back(argv[argn]);
        //program mode command-line switches
        if (find(args.begin(), args.end(), "-?") != args.end()
            || find(args.begin(), args.end(), "-h") != args.end()
            || find(args.begin(), args.end(), "--help") != args.end()) {
            help();
            return 0;
        }

        if (find(args.begin(), args.end(), "--version") != args.end()) {
            version();
            return 0;
        }

        if (find(args.begin(), args.end(), "--fake") != args.end()) {
            generate_fake_data = 1;
            arguments_reflected++;
        }

        if ((op = find(args.begin(), args.end(), "--config-file")) != args.end()
            || (op = find(args.begin(), args.end(), "-f")) != args.end()) {
            op++;
            if (op != args.end()) {
                cfilenameptr = op->c_str();
                arguments_reflected += 2;
                cfile = fopen(cfilenameptr, "r");
                if (!cfile) {
                    fprintf(stderr, "File %s cannot be opened for read.\n", cfilenameptr);
                    return 10;
                } else
                    config_file_given = 1;
            }
        }

        if (find(args.begin(), args.end(), "--write-config") != args.end()
            || find(args.begin(), args.end(), "-w") != args.end()) {
            do_nothing = 1;
            arguments_reflected++;
        }

        if (find(args.begin(), args.end(), "--quiet") != args.end()
            || find(args.begin(), args.end(), "-q") != args.end()
            || find(args.begin(), args.end(), "-q+") != args.end()) {
            be_quiet = 1;
            arguments_reflected++;
        } else {
            if (find(args.begin(), args.end(), "--noquiet") != args.end()
                || find(args.begin(), args.end(), "-q-") != args.end()) {
                be_not_quiet = 1;
                be_quiet = 0;
                arguments_reflected++;
            }
        }
        if ((op = find(args.begin(), args.end(), "-i-")) != args.end()
            || (op = find(args.begin(), args.end(), "--input-")) != args.end()
            || (op = find(args.begin(), args.end(), "--input-file-")) != args.end())
            arguments_reflected++;
        if ((op = find(args.begin(), args.end(), "-i")) != args.end()
            || (op = find(args.begin(), args.end(), "--input")) != args.end()
            || (op = find(args.begin(), args.end(), "--input-file")) != args.end()) {
            op++;
            if (op != args.end()) {
                if (*op != "-") {
                    ifstream *fastatestptr = new ifstream(op->c_str());
                    if (!((*fastatestptr) && fastatestptr->good())) {
                        fprintf(stderr, "Input file %s cannot be opened for read.\n", op->c_str());
                        return 10;
                    } else {
                        fastaptr = fastatestptr;
                        InputFileName = *op;
                    }
                }
                arguments_reflected += 2;
            }
        }

        if ((op = find(args.begin(), args.end(), "-in")) != args.end()
            || (op = find(args.begin(), args.end(), "--input-name")) != args.end()
            || (op = find(args.begin(), args.end(), "--input-file-name")) != args.end()) {
            op++;
            if (op != args.end()) {
                input_file_name_rewritten = 1;
                FakeInputFileName = *op;
                arguments_reflected += 2;
            }
        }
        //switches regulating output
        if ((op = find(args.begin(), args.end(), "--output-file")) != args.end()) {
            op++;
            if (op != args.end()) {
                if (*op != "-") {
                    OutputFileName = *op;
                    ofstream *out_est_ptr = new ofstream(op->c_str(), ios::app);
                    if (!((*out_est_ptr) && out_est_ptr->good())) {
                        fprintf(stderr, "Output file %s cannot be opened for write.\n", op->c_str());
                        return 11;
                    } else {
                        outf_ptr = out_est_ptr;
                        outf_ptr->close();
                        out_ptr = dynamic_cast<ostream *>(outf_ptr);
                        OutputFileOpenedOK = 1;
                    }
                }
                arguments_reflected += 2;
            }
        }

        if ((op = find(args.begin(), args.end(), "--output-fasta")) != args.end()) {
            op++;
            if (op != args.end()) {
                if (*op != "-") {
                    OutputFastAFileName = *op;
                    ofstream *out_est_ptr = new ofstream(op->c_str(), ios::app);
                    if (!((*out_est_ptr) && out_est_ptr->good())) {
                        fprintf(stderr, "Output file %s cannot be opened for write.\n", op->c_str());
                        return 11;
                    } else {
                        outfasta_ptr = out_est_ptr;
                        outfasta_ptr->close();
                        OutputFastAFileOpenedOK = 1;
                    }
                }
                arguments_reflected += 2;
            }
        }

        if ((op = find(args.begin(), args.end(), "--masked-part")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                masked_part_on_output = strtod(value, &ptr);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Masked part value is not double.\n";
                else
                    arguments_reflected += 2;
            }
        }

        if ((op = find(args.begin(), args.end(), "--flanks-length")) != args.end()
            || (op = find(args.begin(), args.end(), "--flanks")) != args.end()
            || (op = find(args.begin(), args.end(), "-k")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                flanks_length = strtoul(value, &ptr, 0);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Flanks length is not an unsigned int.\n";
                else
                    arguments_reflected += 2;
            }
        }

        if ((op = find(args.begin(), args.end(), "--show-spaced-bases")) != args.end()
            || (op = find(args.begin(), args.end(), "--ssb")) != args.end()
            || (op = find(args.begin(), args.end(), "-ssb")) != args.end()) {
            show_letters_in_gap = 1;
            arguments_reflected++;
        } else {
            if ((op = find(args.begin(), args.end(), "--hide-spaced-bases")) != args.end()
                || (op = find(args.begin(), args.end(), "--not-show-spaced-bases")) != args.end()
                || (op = find(args.begin(), args.end(), "--hsb")) != args.end()
                || (op = find(args.begin(), args.end(), "-hsb")) != args.end()) {
                show_letters_in_gap = 0;
                arguments_reflected++;
            }
        }

        if (find(args.begin(), args.end(), "--cgi") != args.end()
            || find(args.begin(), args.end(), "-cgi") != args.end()) {
            output_cgi = 1;
            arguments_reflected++;
        }

        if (find(args.begin(), args.end(), "--xml") != args.end()
            || find(args.begin(), args.end(), "-xml") != args.end()) {
            output_xml = 1;
            arguments_reflected++;
        }

        if (find(args.begin(), args.end(), "--html") != args.end()
            || find(args.begin(), args.end(), "-html") != args.end()) {
            output_html = 1;
            arguments_reflected++;
        }

        if (find(args.begin(), args.end(), "--output-unreliable-result") != args.end()
            || find(args.begin(), args.end(), "--unrel") != args.end()
            || find(args.begin(), args.end(), "-unrel") != args.end()) {
            output_unreliable_result = 1;
            arguments_reflected++;
        }

        if ((op = find(args.begin(), args.end(), "-logfile")) != args.end()
            || (op = find(args.begin(), args.end(), "--logfile")) != args.end()
            || (op = find(args.begin(), args.end(), "-log")) != args.end()
            || (op = find(args.begin(), args.end(), "--log")) != args.end()

        ) {
            op++;
            if (op != args.end()) {
                logfile = *op;
                arguments_reflected += 2;
            }
        }

        if ((op = find(args.begin(), args.end(), "-ip")) != args.end()
            || (op = find(args.begin(), args.end(), "--ip")) != args.end()) {
            op++;
            if (op != args.end()) {
                ip = *op;
                arguments_reflected += 2;
            }
        }

        if ((op = find(args.begin(), args.end(), "-id")) != args.end()
            || (op = find(args.begin(), args.end(), "--id")) != args.end()) {
            op++;
            if (op != args.end()) {
                id = *op;
                arguments_reflected += 2;
            }
        }

        if ((op = find(args.begin(), args.end(), "--time-limit")) != args.end()
            || (op = find(args.begin(), args.end(), "--time")) != args.end()
            || (op = find(args.begin(), args.end(), "-tl")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                time_limit = strtoul(value, &ptr, 0);
                time_limit_is_default = 0;
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Time limit is not an unsigned int.\n";
                else
                    arguments_reflected += 2;
            }
        }

        if (find(args.begin(), args.end(), "--output-tables") != args.end()
            || find(args.begin(), args.end(), "--output-weight-tables") != args.end()
            || find(args.begin(), args.end(), "--pwm") != args.end()
            || find(args.begin(), args.end(), "-pwm") != args.end()) {
            output_weight_tables = 1;
            arguments_reflected++;
        } else {
            if (find(args.begin(), args.end(), "--not-output-tables") != args.end()
                || find(args.begin(), args.end(), "--not-output-weight-tables") != args.end()
                || find(args.begin(), args.end(), "--not-pwm") != args.end()
                || find(args.begin(), args.end(), "--no-pwm") != args.end()
                || find(args.begin(), args.end(), "-not-pwm") != args.end()
                || find(args.begin(), args.end(), "-no-pwm") != args.end()
                || find(args.begin(), args.end(), "--pwm-") != args.end()
                || find(args.begin(), args.end(), "-pwm-") != args.end()) {
                not_output_weight_tables = 1;
                arguments_reflected++;
            }
        }
        if (find(args.begin(), args.end(), "--pwm-ibis-2024") != args.end()
            || find(args.begin(), args.end(), "--pwm-ibis") != args.end()
            || find(args.begin(), args.end(), "--ibis") != args.end()
            || find(args.begin(), args.end(), "-ibis") != args.end()) {
            ibis_output = 1;
            arguments_reflected++;
        }
        if ((op = find(args.begin(), args.end(), "--tf-name")) != args.end()
            || (op = find(args.begin(), args.end(), "--tf")) != args.end()
            || (op = find(args.begin(), args.end(), "-tf-name")) != args.end()
            || (op = find(args.begin(), args.end(), "-tf")) != args.end()
            || (op = find(args.begin(), args.end(), "--TF")) != args.end()
            || (op = find(args.begin(), args.end(), "-TF")) != args.end()) {
            op++;
            if (op != args.end()) {
                TF_name = *op;
                arguments_reflected += 2;
            }
        }
        if ((op = find(args.begin(), args.end(), "--motif-name-suffix")) != args.end()) {
            op++;
            if (op != args.end()) {
                motif_suffix = *op;
                arguments_reflected += 2;
            }
        }
        //
        //mode swithes parsed
        //
        //parsing config file
        int res;
        if (config_file_given) {
            //simple
            if ((res = ReadConfigString(cfile, "mode", mode_string, 99, 0)) < 0)
                return res;
            if (res == 1) {
                config.mode = unknown_mode;
                if (!strncmp(mode_string, "one_strand", 99))
                    config.mode = one_thread;
                if (!strncmp(mode_string, "two_strands", 99))
                    config.mode = two_threads;
                if (config.mode == unknown_mode) {
                    *log_file_ptr << "Unknown mode. I will keep two threads mode.\n";
                    config.mode = two_threads;
                }
            }
            if ((res = ReadConfigString(cfile, "symmetry", mode_string, 99, 0)) < 0)
                return res;
            if (res == 1) {
                config.smode = unknown_smode;
                if (!strncmp(mode_string, "no", 99))
                    config.smode = no;
                if (!strncmp(mode_string, "repeats", 99))
                    config.smode = repeats;
                if (!strncmp(mode_string, "palindromes", 99))
                    config.smode = palindromes;
                if (config.smode == unknown_smode) {
                    *log_file_ptr << "Unknown symmetry. I will keep generic mode.\n";
                    config.smode = no;
                }
            }
            if ((res = ReadConfigDouble(cfile, "motif_absence_prior", &config.motif_absence_prior, 0))
                < 0)
                return res;

            if ((res
                 = ReadConfigBoolean(cfile, "adjust_motif_length", &config.adjust_motif_length, 0))
                < 0)
                return res;

            if ((res = ReadConfigBoolean(cfile, "spaced_motif", &config.spaced_motif, 0)) < 0)
                return res;

            if ((res = ReadConfigString(cfile, "input_caps", mode_string, 99, 0)) < 0)
                return res;

            if (res == 1) {
                if (!strncmp(mode_string, "no", 99) || !strncmp(mode_string, "off", 99)
                    || !strncmp(mode_string, "neglect", 99)) {
                    config.initial_caps = off;
                    config.annealing_caps = off;
                    config.refinement_caps = off;
                } else if (!strncmp(mode_string, "initial", 99)) {
                    config.initial_caps = one;
                    config.annealing_caps = off;
                    config.refinement_caps = off;
                } else if (!strncmp(mode_string, "annealing", 99)) {
                    config.initial_caps = one;
                    config.annealing_caps = one;
                    config.refinement_caps = off;
                } else if (!strncmp(mode_string, "obligatory", 99)) {
                    config.initial_caps = one;
                    config.annealing_caps = one;
                    config.refinement_caps = one;
                } else {
                    *log_file_ptr
                        << "Unknown (deprecated tag input_caps) caps mode. I will keep off.\n";
                    config.initial_caps = off;
                    config.annealing_caps = off;
                    config.refinement_caps = off;
                }
            }

            if ((res = ReadConfigString(cfile, "caps", mode_string, 99, 0)) < 0)
                return res;

            if (res == 1) {
                int decipher_res = config.decipher_caps_mode(mode_string);
                if (decipher_res == 1)
                    *log_file_ptr << "Unknown caps mode. No change.\n";
                if (decipher_res == 2)
                    *log_file_ptr << "Incorrect caps mode sequence. No change.\n";
            }

            if ((res = ReadConfigBoolean(cfile,
                                         "adaptive_pseudocouts",
                                         &config.adaptive_pseudocouts,
                                         0))
                < 0)
                return res;

            if ((res = ReadConfigBoolean(cfile, "common_background", &config.common_background, 0))
                < 0)
                return res;

            if ((res = ReadConfigDouble(cfile, "background_A", &background_A, 0)) < 0)
                return res;
            if (res && !background_A_set) {
                given_background_counters++;
                background_A_set = 1;
            }

            if ((res = ReadConfigDouble(cfile, "background_T", &background_T, 0)) < 0)
                return res;
            if (res && !background_T_set) {
                given_background_counters++;
                background_T_set = 1;
            }

            if ((res = ReadConfigDouble(cfile, "background_G", &background_G, 0)) < 0)
                return res;
            if (res && !background_G_set) {
                given_background_counters++;
                background_G_set = 1;
            }

            if ((res = ReadConfigDouble(cfile, "background_C", &background_C, 0)) < 0)
                return res;
            if (res && !background_C_set) {
                given_background_counters++;
                background_C_set = 1;
            }

            if ((res = ReadConfigUnsignedInt(cfile, "motif_length", &config.motif_length, 0)) < 0)
                return res;

            if ((res = ReadConfigUnsignedInt(cfile,
                                             "maximal_motif_length",
                                             &config.maximal_motif_length,
                                             0))
                < 0)
                return res;

            if ((res = ReadConfigUnsignedInt(cfile,
                                             "minimal_motif_length",
                                             &config.minimal_motif_length,
                                             0))
                < 0)
                return res;

            if ((res = ReadConfigString(cfile, "retrieve_other_sites", mode_string, 99, 0)) < 0)
                return res;

            if (res == 1) {
                if (!strncmp(mode_string, "no", 99) || !strncmp(mode_string, "off", 99))
                    config.retrieve_mode = no_retrieve;
                else if (!strncmp(mode_string, "yes", 99) || !strncmp(mode_string, "on", 99))
                    config.retrieve_mode = yes_retrieve;
                else if (!strncmp(mode_string, "smart", 99))
                    config.retrieve_mode = smart_retrieve;
                else if (!strncmp(mode_string, "pure_smart", 99))
                    config.retrieve_mode = pure_smart_retrieve;
                else {
                    *log_file_ptr << "Unknown post-retrieving instruction. I will keep it off.\n";
                    config.retrieve_mode = no_retrieve;
                }
            }

            if ((res = ReadConfigUnsignedInt(cfile,
                                             "retrieve_other_sites_threshold",
                                             &config.retrieve_threshold,
                                             0))
                < 0)
                return res;

            if ((res = ReadConfigBoolean(cfile, "trim_edges", &config.trim_edges, 0)) < 0)
                return res;

            if ((res = ReadConfigBoolean(cfile, "dumb_trim_edges", &config.trim_edges, 0)) < 0)
                return res;

            //advanced
            if ((res = ReadConfigUnsignedLong(cfile,
                                              "local_step_cycles_between_adjustments",
                                              &config.local_step_cycles_between_adjustments,
                                              0))
                < 0)
                return res;

            if ((res = ReadConfigBoolean(cfile,
                                         "adjustments_during_annealing",
                                         &config.adjustments_during_annealing,
                                         0))
                < 0)
                return res;

            if ((res = ReadConfigString(cfile, "annealing_criterion", mode_string, 99, 0)) < 0)
                return res;
            if (res == 1) {
                config.ctmode = unknown_ctmode;
                if (!strncmp(mode_string, "site_positions", 99))
                    config.ctmode = site_positions;
                if (!strncmp(mode_string, "information", 99))
                    config.ctmode = information;
                if (!strncmp(mode_string, "correlation", 99))
                    config.ctmode = correlation;
                if (config.ctmode == unknown_ctmode) {
                    *log_file_ptr
                        << "Unknown annealing criterion. I will keep the information-based one.\n";
                    config.ctmode = information;
                }
            }
            if ((res = ReadConfigUnsignedInt(cfile,
                                             "cycles_with_minor_change_to_say_cold",
                                             &config.cycles_with_minor_change_to_say_cold,
                                             0))
                < 0)
                return res;

            if ((res = ReadConfigDouble(cfile, "minor_change_level", &config.minor_change_level, 0))
                < 0)
                return res;

            if ((res = ReadConfigUnsignedInt(cfile,
                                             "annealings_number_maximum_is_to_be_global_for",
                                             &config.annealings_number_maximum_is_to_be_global_for,
                                             0))
                < 0)
                return res;

            if ((res = ReadConfigUnsignedLong(cfile,
                                              "steps_number_maximum_is_to_be_global_for",
                                              &config.steps_number_maximum_is_to_be_global_for,
                                              0))
                < 0)
                return res;

            if ((res = ReadConfigBoolean(cfile, "slow_optimisation", &config.slow_optimisation, 0))
                < 0)
                return res;

            if ((res
                 = ReadConfigUnsignedInt(cfile, "chain_fails_after", &config.chain_fails_after, 0))
                < 0)
                return res;

            if ((res = ReadConfigUnsignedInt(cfile, "chains_to_try", &config.chains_to_try, 0)) < 0)
                return res;

            if ((res = ReadConfigUnsignedLong(cfile,
                                              "cycles_per_annealing_attempt",
                                              &config.cycles_per_annealing_attempt,
                                              0))
                < 0)
                return res;

            if ((res = ReadConfigUnsignedInt(cfile,
                                             "annealing_attempts",
                                             &config.annealing_attempts,
                                             0))
                < 0)
                return res;

            if ((res = ReadConfigUnsignedLong(cfile, "random_seed_1", &config.seed1, 0)) < 0)
                return res;

            if ((res = ReadConfigUnsignedLong(cfile, "random_seed_2", &config.seed2, 0)) < 0)
                return res;
        }
        //config file is parsed
        //simple command-line switches

        if (find(args.begin(), args.end(), "--one-strand") != args.end()
            || find(args.begin(), args.end(), "-1") != args.end()) {
            config.mode = one_thread;
            arguments_reflected++;
        } else //we avoid multiple mode definition : olnly one is nparsed
        {
            if (find(args.begin(), args.end(), "--two-strands") != args.end()
                || find(args.begin(), args.end(), "-2") != args.end()) {
                config.mode = two_threads;
                arguments_reflected++;
            }
        }

        if (find(args.begin(), args.end(), "--no-symmetry") != args.end()) {
            config.smode = no;
            arguments_reflected++;
        } else //we avoid multiple mode definition : olnly one is nparsed
        {
            if (find(args.begin(), args.end(), "-tr") != args.end()
                || find(args.begin(), args.end(), "-rp") != args.end()
                || find(args.begin(), args.end(), "--tandems") != args.end()
                || find(args.begin(), args.end(), "--repeats") != args.end()) {
                config.smode = repeats;
                arguments_reflected++;
            } else {
                if (find(args.begin(), args.end(), "--palindromes") != args.end()
                    || find(args.begin(), args.end(), "-o") != args.end()) {
                    config.smode = palindromes;
                    arguments_reflected++;
                }
            }
        }

        if ((op = find(args.begin(), args.end(), "--absence-prior")) != args.end()
            || (op = find(args.begin(), args.end(), "-p")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                config.motif_absence_prior = strtod(value, &ptr);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Absence prior value is not double.\n";
                else
                    arguments_reflected += 2;
            }
        }

        if (find(args.begin(), args.end(), "--not-adjust-length") != args.end()
            || find(args.begin(), args.end(), "-a-") != args.end()) {
            config.adjust_motif_length = 0;
            arguments_reflected++;
        } else //we avoid -a+ -a- sutation - one of the is unparsed
        {
            if (find(args.begin(), args.end(), "--adjust-length") != args.end()
                || find(args.begin(), args.end(), "-a+") != args.end()) {
                config.adjust_motif_length = 1;
                arguments_reflected++;
            }
        }

        if (find(args.begin(), args.end(), "--no-spacers") != args.end()
            || find(args.begin(), args.end(), "--not-spaced") != args.end()) {
            config.spaced_motif = 0;
            arguments_reflected++;
        } else //we avoid -a+ -a- sutation - one of the is unparsed
        {
            if (find(args.begin(), args.end(), "--spacers") != args.end()
                || find(args.begin(), args.end(), "--spaced") != args.end()) {
                config.spaced_motif = 1;
                arguments_reflected++;
            }
        }

        if ((op = find(args.begin(), args.end(), "--length")) != args.end()
            || (op = find(args.begin(), args.end(), "-l")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                config.motif_length = strtoul(value, &ptr, 0);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Motif length value is not unsigned int.\n";
                else
                    arguments_reflected += 2;
            }
        }

        if ((op = find(args.begin(), args.end(), "--A")) != args.end()
            || (op = find(args.begin(), args.end(), "--background-A")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                background_A = strtod(value, &ptr);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Background_A is not double.\n";
                else {
                    if (!background_A_set) {
                        background_A_set = 1;
                        given_background_counters++;
                    }
                    arguments_reflected += 2;
                }
            }
        }

        if ((op = find(args.begin(), args.end(), "--T")) != args.end()
            || (op = find(args.begin(), args.end(), "--background-T")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                background_T = strtod(value, &ptr);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Background_T is not double.\n";
                else {
                    if (!background_T_set) {
                        background_T_set = 1;
                        given_background_counters++;
                    }
                    arguments_reflected += 2;
                }
            }
        }

        if ((op = find(args.begin(), args.end(), "--G")) != args.end()
            || (op = find(args.begin(), args.end(), "--background-G")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                background_G = strtod(value, &ptr);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Background_G is not double.\n";
                else {
                    if (!background_G_set) {
                        background_G_set = 1;
                        given_background_counters++;
                    }
                    arguments_reflected += 2;
                }
            }
        }

        if ((op = find(args.begin(), args.end(), "--C")) != args.end()
            || (op = find(args.begin(), args.end(), "--background-C")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                background_C = strtod(value, &ptr);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Background_C is not double.\n";
                else {
                    if (!background_C_set) {
                        background_C_set = 1;
                        given_background_counters++;
                    }
                    arguments_reflected += 2;
                }
            }
        }

        if (find(args.begin(), args.end(), "--footprint") != args.end()
            || find(args.begin(), args.end(), "--footprints") != args.end()
            || find(args.begin(), args.end(), "-fp") != args.end()
            || find(args.begin(), args.end(), "--fp") != args.end()) {
            config.common_background = 1;
            config.refinement_caps = one;
            config.annealing_caps = one;
            config.initial_caps = one;
            arguments_reflected++;
        }

        if (find(args.begin(), args.end(), "--footprinta") != args.end()
            || find(args.begin(), args.end(), "--footprintsa") != args.end()
            || find(args.begin(), args.end(), "-fpa") != args.end()
            || find(args.begin(), args.end(), "--fpa") != args.end()) {
            config.common_background = 1;
            config.refinement_caps = one;
            config.annealing_caps = one;
            config.initial_caps = all;
            arguments_reflected++;
        }

        if ((op = find(args.begin(), args.end(), "--caps")) != args.end()
            || (op = find(args.begin(), args.end(), "-caps")) != args.end()
            || (op = find(args.begin(), args.end(), "--caps-mode")) != args.end()
            || (op = find(args.begin(), args.end(), "-caps-mode")) != args.end()) {
            op++;
            if (op != args.end()) {
                int decipher_res = config.decipher_caps_mode(*op);
                if (decipher_res == 1)
                    *log_file_ptr << "Unknown caps mode. No change.\n";
                if (decipher_res == 2)
                    *log_file_ptr << "Incorrect caps mode sequence. No change.\n";
                if (decipher_res == 0)
                    arguments_reflected += 2;
            }
        } else if (find(args.begin(), args.end(), "--initial-input-caps") != args.end()
                   || find(args.begin(), args.end(), "-iics") != args.end()
                   || find(args.begin(), args.end(), "--iics") != args.end()) {
            config.initial_caps = one;
            config.annealing_caps = off;
            config.refinement_caps = off;
            arguments_reflected++;
        } else if (find(args.begin(), args.end(), "--annealing-input-caps") != args.end()
                   || find(args.begin(), args.end(), "-aics") != args.end()
                   || find(args.begin(), args.end(), "--aics") != args.end()) {
            config.initial_caps = one;
            config.annealing_caps = one;
            config.refinement_caps = off;
            arguments_reflected++;
        } else if (find(args.begin(), args.end(), "--obligatory-input-caps") != args.end()
                   || find(args.begin(), args.end(), "-oics") != args.end()
                   || find(args.begin(), args.end(), "--oics") != args.end()) {
            config.initial_caps = one;
            config.annealing_caps = one;
            config.refinement_caps = one;
            arguments_reflected++;
        } else if (find(args.begin(), args.end(), "--neglect-input-caps") != args.end()
                   || find(args.begin(), args.end(), "-nics") != args.end()
                   || find(args.begin(), args.end(), "--nics") != args.end()) {
            config.initial_caps = off;
            config.annealing_caps = off;
            config.refinement_caps = off;
            arguments_reflected++;
        };

        if ((find(args.begin(), args.end(), "--adaptive-pseudocouts") != args.end())
            || (find(args.begin(), args.end(), "-adap") != args.end())
            || (find(args.begin(), args.end(), "-adaptive-pseudocouts") != args.end())
            || (find(args.begin(), args.end(), "--adap") != args.end())) {
            config.adaptive_pseudocouts = 1;
            arguments_reflected++;
        } else //we avoid -a+ -a- sutation - one of the is unparsed
        {
            if ((find(args.begin(), args.end(), "--no-adaptive-pseudocouts") != args.end())
                || (find(args.begin(), args.end(), "--nadap") != args.end())
                || (find(args.begin(), args.end(), "-no-adaptive-pseudocouts") != args.end())
                || (find(args.begin(), args.end(), "-nadap") != args.end())) {
                config.adaptive_pseudocouts = 0;
                arguments_reflected++;
            }
        }

        if (find(args.begin(), args.end(), "--common-background") != args.end()) {
            config.common_background = 1;
            arguments_reflected++;
        } else //we avoid -a+ -a- sutation - one of the is unparsed
        {
            if (find(args.begin(), args.end(), "--no-common-background") != args.end()) {
                config.common_background = 0;
                arguments_reflected++;
            }
        }

        if ((op = find(args.begin(), args.end(), "--maximal-length")) != args.end()
            || (op = find(args.begin(), args.end(), "--max-length")) != args.end()
            || (op = find(args.begin(), args.end(), "-l+")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                config.maximal_motif_length = strtoul(value, &ptr, 0);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Maximal motif length value is not unsigned int.\n";
                else
                    arguments_reflected += 2;
            }
        }

        if ((op = find(args.begin(), args.end(), "--minimal-length")) != args.end()
            || (op = find(args.begin(), args.end(), "--min-length")) != args.end()
            || (op = find(args.begin(), args.end(), "-l-")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                config.minimal_motif_length = strtoul(value, &ptr, 0);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Minimal motif length value is not unsigned int.\n";
                else
                    arguments_reflected += 2;
            }
        }

        if (find(args.begin(), args.end(), "--not-retrieve-other-sites") != args.end()
            || find(args.begin(), args.end(), "-no-retrieve") != args.end()
            || find(args.begin(), args.end(), "--no-retrieve") != args.end()
            || find(args.begin(), args.end(), "-r-") != args.end()) {
            config.retrieve_mode = no_retrieve;
            arguments_reflected++;
        } else if (find(args.begin(), args.end(), "--smart-retrieve-other-sites") != args.end()
                   || find(args.begin(), args.end(), "--retrieve-smart") != args.end()
                   || find(args.begin(), args.end(), "-rs") != args.end()) {
            config.retrieve_mode = smart_retrieve;
            arguments_reflected++;
        } else if (find(args.begin(), args.end(), "--pure-smart-retrieve-other-sites") != args.end()
                   || find(args.begin(), args.end(), "--retrieve-pure-smart") != args.end()
                   || find(args.begin(), args.end(), "-rps") != args.end()) {
            config.retrieve_mode = pure_smart_retrieve;
            arguments_reflected++;
        }

        else //we avoid -r+ -r- sutation - one of the is unparsed
        {
            if (find(args.begin(), args.end(), "--retrieve-other-sites") != args.end()
                || find(args.begin(), args.end(), "-r+") != args.end()) {
                config.retrieve_mode = yes_retrieve;
                arguments_reflected++;
            }
        }

        if ((op = find(args.begin(), args.end(), "--retrieve-other-sites-threshold")) != args.end()
            || (op = find(args.begin(), args.end(), "-r")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                config.retrieve_threshold = strtoul(value, &ptr, 0);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Number of the threshold site is not an unsigned int.\n";
                else
                    arguments_reflected += 2;
            }
        }

        if (find(args.begin(), args.end(), "--not-trim-edges") != args.end()
            || find(args.begin(), args.end(), "--not-trim") != args.end()
            || find(args.begin(), args.end(), "-not-trim") != args.end()
            || find(args.begin(), args.end(), "--not-trim-edges") != args.end()
            || find(args.begin(), args.end(), "--trim-") != args.end()
            || find(args.begin(), args.end(), "-trim-") != args.end()) {
            config.trim_edges = 0;
            arguments_reflected++;
        } else //we avoid --trim --not-trim sutation - one of the is unparsed
        {
            if (find(args.begin(), args.end(), "--trim-edges") != args.end()
                || find(args.begin(), args.end(), "--trim") != args.end()
                || find(args.begin(), args.end(), "-trim-edges") != args.end()
                || find(args.begin(), args.end(), "-trim") != args.end()) {
                config.trim_edges = 1;
                arguments_reflected++;
            }
        }

        if (find(args.begin(), args.end(), "--not-dumb-trim-edges") != args.end()
            || find(args.begin(), args.end(), "--not-dumb-trim") != args.end()
            || find(args.begin(), args.end(), "-not-dumb-trim") != args.end()
            || find(args.begin(), args.end(), "--not-dumb-trim-edges") != args.end()
            || find(args.begin(), args.end(), "--dumb-trim-") != args.end()
            || find(args.begin(), args.end(), "-dumb-trim-") != args.end()) {
            config.dumb_trim_edges = 0;
            arguments_reflected++;
        } else //we avoid --dumb-trim --not-dumb-trim sutation - one of the is unparsed
        {
            if (find(args.begin(), args.end(), "--dumb-trim-edges") != args.end()
                || find(args.begin(), args.end(), "--dumb-trim") != args.end()
                || find(args.begin(), args.end(), "-dumb-trim-edges") != args.end()
                || find(args.begin(), args.end(), "-dumb-trim") != args.end()) {
                config.dumb_trim_edges = 1;
                arguments_reflected++;
            }
        }

        //advanced command-line switches
        if ((op = find(args.begin(), args.end(), "--local-cycles")) != args.end()
            || (op = find(args.begin(), args.end(), "-e")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                config.local_step_cycles_between_adjustments = strtoul(value, &ptr, 0);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr
                        << "Local cycles between adjustments value is not unsigned long.\n";
                else
                    arguments_reflected += 2;
            }
        }

        if (find(args.begin(), args.end(), "--annealing-adjustments") != args.end()
            || find(args.begin(), args.end(), "-n+") != args.end()) {
            config.adjustments_during_annealing = 1;
            arguments_reflected++;
        } else //we avoid -n+ -n- sutation - one of the is unparsed
        {
            if (find(args.begin(), args.end(), "--not-annealing-adjustments") != args.end()
                || find(args.begin(), args.end(), "-n-") != args.end()) {
                config.adjustments_during_annealing = 0;
                arguments_reflected++;
            }
        }

        if (find(args.begin(), args.end(), "--site-positions-annealing-criterion") != args.end()
            || find(args.begin(), args.end(), "-sac") != args.end()) {
            config.ctmode = site_positions;
            arguments_reflected++;
        } else //we avoid multiple mode definition : olnly one is nparsed
        {
            if (find(args.begin(), args.end(), "--information-annealing-criterion") != args.end()
                || find(args.begin(), args.end(), "-iac") != args.end()) {
                config.ctmode = information;
                arguments_reflected++;
            } else {
                if (find(args.begin(), args.end(), "--correlation-annealing-criterion") != args.end()
                    || find(args.begin(), args.end(), "-cac") != args.end()) {
                    config.ctmode = correlation;
                    arguments_reflected++;
                }
            }
        }

        if ((op = find(args.begin(), args.end(), "--cycles-to-say-cold")) != args.end()
            || (op = find(args.begin(), args.end(), "-y")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                config.cycles_with_minor_change_to_say_cold = strtoul(value, &ptr, 0);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Local cycles to stop annealing value is not unsigned long.\n";
                else
                    arguments_reflected += 2;
            }
        }

        if ((op = find(args.begin(), args.end(), "--annealing-stop-level")) != args.end()
            || (op = find(args.begin(), args.end(), "-v")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                config.minor_change_level = strtod(value, &ptr);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Annealing stop bound level not double.\n";
                else
                    arguments_reflected += 2;
            }
        }

        if ((op = find(args.begin(), args.end(), "--global-area-annealings")) != args.end()
            || (op = find(args.begin(), args.end(), "-g")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                config.annealings_number_maximum_is_to_be_global_for = strtoul(value, &ptr, 0);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Annealings duration maximum is to hold to be global value is "
                                     "not unsigned int.\n";
                else
                    arguments_reflected += 2;
            }
        }

        if ((op = find(args.begin(), args.end(), "--global-area-steps")) != args.end()
            || (op = find(args.begin(), args.end(), "-t")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                config.steps_number_maximum_is_to_be_global_for = strtoul(value, &ptr, 0);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Steps naumber maximum is to hold to be global value is not "
                                     "unsigned long.\n";
                else
                    arguments_reflected += 2;
            }
        }

        if (find(args.begin(), args.end(), "--slow-optimisation") != args.end()
            || find(args.begin(), args.end(), "--slow-optimization") != args.end()
            || find(args.begin(), args.end(), "-s") != args.end()) {
            config.slow_optimisation = 1;
            arguments_reflected++;
        }

        if ((op = find(args.begin(), args.end(), "--chain-fails-after")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                config.chain_fails_after = strtoul(value, &ptr, 0);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr
                        << "Number of forbidden adjsjmts to fail chain is not unsigned int.\n";
                else
                    arguments_reflected += 2;
            }
        }

        if ((op = find(args.begin(), args.end(), "--chains-to-try")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                config.chains_to_try = strtoul(value, &ptr, 0);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Number of chains to try is not unsigned int.\n";
                else
                    arguments_reflected += 2;
            }
        }

        if ((op = find(args.begin(), args.end(), "--cycles_per_annealing_attempt")) != args.end()
            || (op = find(args.begin(), args.end(), "--ac")) != args.end()
            || (op = find(args.begin(), args.end(), "-ac")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                config.cycles_per_annealing_attempt = strtoul(value, &ptr, 0);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr
                        << "Number of cycles per annealing attempt is not unsigned int.\n";
                else
                    arguments_reflected += 2;
            }
        }

        if ((op = find(args.begin(), args.end(), "--annealing-attempts")) != args.end()
            || (op = find(args.begin(), args.end(), "--at")) != args.end()
            || (op = find(args.begin(), args.end(), "-at")) != args.end()) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                config.annealing_attempts = strtoul(value, &ptr, 0);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Number of annealing attempts is not unsigned int.\n";
                else
                    arguments_reflected += 2;
            }
        }

        if ((op = find(args.begin(), args.end(), "--random-seed-1")) != args.end()
            || (op = find(args.begin(), args.end(), "--seed-1")) != args.end()
            || (op = find(args.begin(), args.end(), "--seed1")) != args.end()
            || (op = find(args.begin(), args.end(), "-seed-1")) != args.end()
            || (op = find(args.begin(), args.end(), "-seed1")) != args.end()

        ) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                config.seed1 = strtoul(value, &ptr, 0);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Seed 1 is not an unsigned long.\n";
                else
                    arguments_reflected += 2;
            }
        }

        if ((op = find(args.begin(), args.end(), "--random-seed-2")) != args.end()
            || (op = find(args.begin(), args.end(), "--seed-2")) != args.end()
            || (op = find(args.begin(), args.end(), "--seed2")) != args.end()
            || (op = find(args.begin(), args.end(), "-seed-2")) != args.end()
            || (op = find(args.begin(), args.end(), "-seed2")) != args.end()

        ) {
            op++;
            if (op != args.end()) {
                value = op->c_str();
                config.seed2 = strtoul(value, &ptr, 0);
                if (ptr != (value + strlen(value)))
                    *log_file_ptr << "Seed 2 is not an unsigned long.\n";
                else
                    arguments_reflected += 2;
            }
        }
        //test woraround
        if ((find(args.begin(), args.end(), "--no-permutations") != args.end())
            || (find(args.begin(), args.end(), "--no-permut") != args.end())) {
            forbid_permutations = 1;
            arguments_reflected++;
        }
        //testing
        if (args.size() != arguments_reflected) {
            *log_file_ptr << "It seemes to me that some arguments are unknown or duplicated..."
                          << endl
                          << "I will output the configuration if you give me data and exit." << endl
                          << endl;
            do_nothing = 1;
            errlevel = 10;
        }

        //trim desision
        if (config.trim_edges)
            config.dumb_trim_edges = 0;
        //making decisions about annealing test

        if (config.ctmode == information) {
            if (config.minor_change_level < 0)
                config.minor_change_level = 0.1;
            if (config.cycles_with_minor_change_to_say_cold == 0)
                config.cycles_with_minor_change_to_say_cold = 5;
        }
        if (config.ctmode == site_positions) {
            if (config.minor_change_level < 0)
                config.minor_change_level = 0.7;
            if (config.cycles_with_minor_change_to_say_cold == 0)
                config.cycles_with_minor_change_to_say_cold = 5;
        }
        //making desisions about timer
        if (time_limit_is_default && config.slow_optimisation)
            time_limit *= 10;
        //making decisions about input/output

        if (ibis_output) {
            //exclusive ibis mode
            if (TF_name == "noname") {
                if (InputFileName == "") {
                    *log_file_ptr << "Ibis and no tf name!\n";
                } else {
                    TF_name = InputFileName;
                    //TF_name=std::filesystem::path(InputFileName).stem();
                }
            }
            config.minimal_motif_length = 5;
            config.maximal_motif_length = 30;
        }

        if (do_nothing)
            be_quiet = 1;

        if (output_cgi) {
            if (!be_not_quiet)
                be_quiet = 1;
            if (!not_output_weight_tables)
                output_weight_tables = 1;
            log_file_ptr = &cout;
            output_html = 1;
        }
        if (output_html || output_xml) {
            if (!be_not_quiet)
                be_quiet = 1;
        }

        if (output_html)
            output_xml = 0;

        if (input_file_name_rewritten)
            InputFileName = FakeInputFileName;

        if (config.retrieve_mode == undef_retrieve) //default
            config.retrieve_mode = (config.slow_optimisation ? no_retrieve : smart_retrieve);

        if (config.retrieve_threshold == 0x10000) //default
        {
            if (config.retrieve_mode == no_retrieve || config.retrieve_mode == pure_smart_retrieve)
                config.retrieve_threshold = 0;
            else
                config.retrieve_threshold = 1;
        }

        if (config.retrieve_mode == no_retrieve and (config.retrieve_threshold != 0)) {
            *log_file_ptr
                << "Non-zero threshold is given for no posprocessing... I put it back to 0.\n";
            config.retrieve_threshold = 0;
        }

        if (config.retrieve_mode == pure_smart_retrieve and (config.retrieve_threshold != 0)) {
            *log_file_ptr
                << "Non-zero threshold is given for \"pure smart\" ... I put it back to 0.\n";
            config.retrieve_threshold = 0;
        }

        if (config.retrieve_mode == smart_retrieve and (config.retrieve_threshold == 0)) {
            *log_file_ptr << "Zero threshold for \"smart\" ... It is \"pure smart\", actually.\n";
        }

        //decisions about input/output are made

        //deciding what to do with the background is we have at least
        //one given background counter, otherwise let background() to
        //be empty and the sampler will not use it.
        if (given_background_counters) {
            if (given_background_counters != 4) {
                *log_file_ptr << "You have given only " << given_background_counters
                              << " background counters." << endl
                              << endl;
                do_nothing = 1;
                errlevel = 10;
            } else //4
            {
                config.background.resize(4);
                config.background[0] = background_A;
                config.background[1] = background_T;
                config.background[2] = background_G;
                config.background[3] = background_C;
                if (find_if(config.background.begin(),
                            config.background.end(),
                            bind(less_equal<double>(), placeholders::_1, 0))
                    != config.background.end()) {
                    *log_file_ptr << "You have given me a nonpositive background counter." << endl
                                  << endl;
                    do_nothing = 1;
                    errlevel = 10;
                } else
                    config.common_background = 1;
            }
        }

        //not to gett too much stupid output
        if (do_nothing)
            be_quiet = 1;
        //init random generator

        //
        rinit(config.seed1, config.seed2);

        if (!be_quiet)
            *log_file_ptr << "Random generator is created...\n" << flush;

        //random generator is started
        SequencesPile sp;

        if (generate_fake_data) {
            if (!config.spaced_motif)
                createFakeExperimentalBunch(sp, be_quiet, log_file_ptr);
            else
                createFakeExperimentalGappedBunch(sp, be_quiet, log_file_ptr);
            InputFileName = "fake generated data";
            //fake *sb is ready
        } else {
            if (!be_quiet)
                *log_file_ptr << "reading FastA\n";
            try {
                *fastaptr >> sp;
            } catch (const DumbException &de) {
                *log_file_ptr << de;
                return 111;
            }
            *log_file_ptr << sp.Read_Diagnostics;
            // *sb is ready
        }

        if (sp.size() > 1) {
            if (!be_quiet)
                *log_file_ptr << "Sequence bunch is read...\n" << flush;
        } else {
            *log_file_ptr << "One sequence or nothing is given...\n" << flush;
            return 17;
        }

        if (sp.size() <= 1) {
            *log_file_ptr << "One good sequence or no good sequences is given...\n" << flush;
            return 17;
        }
        try {
            switch (config.smode) {
            case no:
                switch (config.mode) {
                case one_thread:
                    sampler = new LookingForAtgcMotifsInOneThreadMultinomialGibbs(
                        sp,
                        config.motif_absence_prior,
                        config.local_step_cycles_between_adjustments,
                        config.initial_caps,
                        config.annealing_caps,
                        config.refinement_caps,
                        config.adaptive_pseudocouts,
                        config.common_background,
                        config.background,
                        be_quiet,
                        *log_file_ptr);
                    break;
                case two_threads:
                    sampler = new LookingForAtgcMotifsInTwoThreadsMultinomialGibbs(
                        sp,
                        config.motif_absence_prior,
                        config.local_step_cycles_between_adjustments,
                        config.initial_caps,
                        config.annealing_caps,
                        config.refinement_caps,
                        config.adaptive_pseudocouts,
                        config.common_background,
                        config.background,
                        be_quiet,
                        *log_file_ptr);
                    break;
                case unknown_mode:
                    break;
                }
                break;
            case repeats:
                switch (config.mode) {
                case one_thread:
                    sampler = new LookingForAtgcRepeatsInOneThreadMultinomialGibbs(
                        sp,
                        config.motif_absence_prior,
                        config.local_step_cycles_between_adjustments,
                        config.initial_caps,
                        config.annealing_caps,
                        config.refinement_caps,
                        config.adaptive_pseudocouts,
                        config.common_background,
                        config.background,
                        be_quiet,
                        *log_file_ptr);
                    break;
                case two_threads:
                    sampler = new LookingForAtgcRepeatsInTwoThreadsMultinomialGibbs(
                        sp,
                        config.motif_absence_prior,
                        config.local_step_cycles_between_adjustments,
                        config.initial_caps,
                        config.annealing_caps,
                        config.refinement_caps,
                        config.adaptive_pseudocouts,
                        config.common_background,
                        config.background,
                        be_quiet,
                        *log_file_ptr);
                    break;
                case unknown_mode:
                    break;
                }
                break;
            case palindromes:
                sampler = new LookingForAtgcPalindromesMultinomialGibbs(
                    sp,
                    config.motif_absence_prior,
                    config.local_step_cycles_between_adjustments,
                    config.initial_caps,
                    config.annealing_caps,
                    config.refinement_caps,
                    config.adaptive_pseudocouts,
                    config.common_background,
                    config.background,
                    be_quiet,
                    *log_file_ptr);
                break;
            case unknown_smode:;
            }
        } catch (DumbException &de) {
            *log_file_ptr << de;
            return 10;
        };

        /*retrieving and setting bound values from and fro the sampler*/

        sampler->SayLengthBounds(config.minimal_motif_length,
                                 config.motif_length,
                                 config.maximal_motif_length);

        if (config.minimal_motif_length > config.maximal_motif_length) {
            *log_file_ptr << "The length requirements from the data contradict with the\n" << flush;
            *log_file_ptr << "the length parameters.\n" << flush;
            return 31;
        }

        if (config.minimal_motif_length == config.maximal_motif_length
            && config.adjust_motif_length) {
            config.adjust_motif_length = 0;
        }

        config.motif_absence_prior = sampler->motif_absence_prior;

        /*****************************************************************************/
        /* doing nothing (--write-config)                                            */
        /*****************************************************************************/
        if (do_nothing) {
            /* retrieving the parameters values from the sampler                         */
            config.print(cout);
            return errlevel;
        }
        /*****************************************************************************/

        //test workaround
        if (forbid_permutations)
            sampler->permutation_forbidden = 1;

        LogRecorder logger(logfile, id, ip);
        if (!be_quiet)
            *log_file_ptr << "Sampler is created..." << flush;

        if (config.slow_optimisation && !config.adjust_motif_length) {
            diagnostics
                << "The slow optimisation mode makes no sense when the site length is pre-given.\n";
            config.slow_optimisation = 0;
        }
        if (!config.slow_optimisation) {
            try {
                sampler->find_maximum(config.motif_length,
                                      config.minimal_motif_length,
                                      config.maximal_motif_length,
                                      config.adjust_motif_length,
                                      config.spaced_motif,
                                      config.steps_number_maximum_is_to_be_global_for,
                                      config.annealings_number_maximum_is_to_be_global_for,
                                      config.cycles_with_minor_change_to_say_cold,
                                      config.minor_change_level,
                                      config.ctmode,
                                      config.adjustments_during_annealing,
                                      config.chain_fails_after,
                                      config.chains_to_try,
                                      config.cycles_per_annealing_attempt,
                                      config.annealing_attempts,
                                      time_limit,
                                      logger);
            } catch (LookingForAtgcMotifsMultinomialGibbs::LostInSpaceOnSecondaryAnnealing &lo) {
                diagnostics.status = fatal;
                diagnostics << lo;
            } catch (LookingForAtgcMotifsMultinomialGibbs::AllInitialAnnealingAttemptsFailed &lo) {
                diagnostics.status = fatal;
                diagnostics << lo;
            } catch (LookingForAtgcMotifsMultinomialGibbs::LostInSpaceOnTracing &lo) {
                diagnostics.status = unreliable;
                diagnostics << lo;
            } catch (LookingForAtgcMotifsMultinomialGibbs::TooMuchChainsFailed &too) {
                diagnostics.status = unreliable;
                diagnostics << too;
            } catch (LookingForAtgcMotifsMultinomialGibbs::TimeLimitException &tle) {
                diagnostics.status = fatal;
                diagnostics << tle;
            } catch (LookingForAtgcMotifsMultinomialGibbs::TooResrtrictiveCapsMode &trcm) {
                diagnostics.status = fatal;
                diagnostics << trcm;
            } catch (DumbException &error) {
                diagnostics.status = fatal;
                diagnostics << error;
            }
        } else {
            try {
                sampler->find_maximum_slowly(config.minimal_motif_length,
                                             config.maximal_motif_length,
                                             config.spaced_motif,
                                             G_values,
                                             config.steps_number_maximum_is_to_be_global_for,
                                             config.annealings_number_maximum_is_to_be_global_for,
                                             config.cycles_with_minor_change_to_say_cold,
                                             config.minor_change_level,
                                             config.ctmode,
                                             config.adjustments_during_annealing,
                                             config.chain_fails_after,
                                             config.chains_to_try,
                                             config.cycles_per_annealing_attempt,
                                             config.annealing_attempts,
                                             time_limit,
                                             logger,
                                             diagnostics);
            } catch (LookingForAtgcMotifsMultinomialGibbs::TimeLimitException &tle) {
                diagnostics.status = fatal;
                diagnostics << tle;
            } catch (LookingForAtgcMotifsMultinomialGibbs::TooResrtrictiveCapsMode &trcm) {
                diagnostics.status = fatal;
                diagnostics << trcm;
            } catch (LookingForAtgcMotifsMultinomialGibbs::SlowSearchFailed &fai) {
                diagnostics.status = fatal;
                diagnostics << fai;
            } catch (DumbException &error) {
                diagnostics.status = fatal;
                diagnostics << error;
            }
        }

        if (sampler->theState.pattern_length() == 0) {
            diagnostics.status = fatal;
            diagnostics << "The state we found has zero site length!\n";
        }

        if (output_html)
            diagnostics.output_mode = html_output;
        else if (output_xml)
            diagnostics.output_mode = xml_output;
        else
            diagnostics.output_mode = txt_output;

        if (diagnostics.status == fatal
            || (diagnostics.status == unreliable && !output_unreliable_result)) //finishing
        {
            if (OutputFileOpenedOK)
                outf_ptr->open(OutputFileName.c_str());
            if (output_html) {
                output_html_header(*out_ptr, InputFileName);
                if (diagnostics.str() != "")
                    *out_ptr << diagnostics;
                config.print(*out_ptr, "<!-- The configuration parameters:\n", "-->\n", "");
                output_html_footer(*out_ptr);
            } else if (output_xml) {
                output_xml_header(*out_ptr);
                if (diagnostics.str() != "")
                    *out_ptr << diagnostics;
                config.print(*out_ptr,
                             "<comment name=\"configuration_parameters\">\n",
                             "</comment>\n",
                             "");
                output_xml_footer(*out_ptr);
            } else {
                output_text_header(*out_ptr, InputFileName);
                if (diagnostics.str() != "")
                    *out_ptr << diagnostics;
                config.print(*out_ptr, "#The configuration parameters:\n", "", "#");
            }
            if (OutputFileOpenedOK)
                outf_ptr->close();
        } else {
            if (!be_quiet)
                *log_file_ptr << "Finished. Giving results.....\n";
            if (!be_quiet)
                *log_file_ptr << sampler->theState;
            vector<Profile *> results; //for future
            //giving the results themself
            Profile profile(sp,
                            sampler->theState,
                            *sampler,
                            config.adjust_motif_length,
                            config.minimal_motif_length,
                            config.maximal_motif_length,
                            config.slow_optimisation,
                            config.retrieve_mode,
                            config.retrieve_threshold,
                            config.trim_edges,
                            config.dumb_trim_edges,
                            config.common_background,
                            config.background,
                            diagnostics,
                            order,
                            flanks_length,
                            show_letters_in_gap,
                            G_values);
            //DEBUG
            //cout<<"\n<P>= "<<diagnostics<<" =<P>\n";
            results.push_back(&profile);
            results[0]->output_additional_information = output_weight_tables;
            try {
                if (OutputFileOpenedOK)
                    outf_ptr->open(OutputFileName.c_str());
                if (output_html) {
                    output_html_header(*out_ptr, InputFileName);
                    results[0]->html_output(*out_ptr);
                    results[0]->short_html_output(*out_ptr);
                    config.print(*out_ptr, "<!-- The configuration parameters:\n", "-->\n", "");
                    output_html_footer(*out_ptr);
                } else if (output_xml) {
                    output_xml_header(*out_ptr);
                    results[0]->xml_output(*out_ptr, id, InputFileName);
                    config.print(*out_ptr,
                                 "<comment name=\"configuration_parameters\">\n",
                                 "</comment>\n",
                                 "");
                    output_xml_footer(*out_ptr);
                } else if (ibis_output) {
                    if (motif_suffix == "noname") {
                        motif_suffix = "_SeSiMCMC_I";
                    }
                    results[0]->ibis_pwm_output(*out_ptr, TF_name, motif_suffix);
                } else //text
                {
                    output_text_header(*out_ptr, InputFileName);
                    results[0]->text_output(*out_ptr);
                    config.print(*out_ptr, "#The configuration parameters:\n", "", "#");
                }
                if (OutputFileOpenedOK)
                    outf_ptr->close();
                if (OutputFastAFileOpenedOK) {
                    outfasta_ptr->open(OutputFastAFileName.c_str());
                    sp.put_mask(*results[0], masked_part_on_output);
                    *outfasta_ptr << sp;
                    outfasta_ptr->close();
                }
            } catch (DumbException &error) {
                *log_file_ptr << "Cannot output report data\n" << error;
                return 1100;
            }
        }

        if (!be_quiet)
            *log_file_ptr << "Deleting sampler...\n";
        delete sampler;
    } //main try-block
    catch (const DumbException &de) {
        std::cerr << de << '\n';
    }

    catch (const std::exception &ex) {
        std::cerr << ex.what();
    }

    catch (...) {
        std::cerr << "Unexpected exception caught\n";
        return -1;
    }
    return 0;
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




void createFakeExperimentalGappedBunch(SequencesPile &sp,
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
	unsigned short tpattern[10]={1,2,3,3,0,0,1,2,3,1};
	unsigned short cpattern[10]={2,4,1,2,0,0,4,4,1,2};
//	unsigned short tpattern[10]={1,2,3,0,0,0,0,2,3,1};
//	unsigned short cpattern[10]={2,4,1,0,0,0,0,4,1,2};
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
				for (unsigned int j=0;j<p_leng;j++)
					if (pattern[j]) a[fake_data_positions[i]+j]=pattern[j];
				  //we do not map zeroes, it is gap
				pattern[ch_pos]=old_let; //restore
			}
			else
			{
				for (unsigned int j=0;j<p_leng;j++)
					if (pattern[j]) a[fake_data_positions[i]+j]=pattern[j];
				  //we do not map zeroes, it is gap
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

