#include <iostream>
#include <iomanip>

#include "../Sequences.hpp"
#include "../SymbolsCounter.hpp"
#include "../MarkovChainState.hpp"

/*
 * 2 3 4 5 6 7 8 9 10 11
 * 1 2 3 4 5 6 7 8  9 10 11 12 13 14 15  ***** current
 * 1 2 3 4 5 6 7 8  9 10 11
 */
int main()
{
	SequencesPile * sp;
	SymbolsCounter * sc;
	MarkovChainState mcs(3); //3 sequences
	unsigned int pat_len=5;
	cout<<"Start"<<endl<<flush;
	unsigned int letters=4;
//	unsigned int current_sequence=2;
	unsigned short *a1,*a2,*a3;
	unsigned int l1=10,l2=15,l3=11,i;
	a1=(unsigned short*)calloc(l1+1,sizeof(unsigned short));
	a2=(unsigned short*)calloc(l2+1,sizeof(unsigned short));
	a3=(unsigned short*)calloc(l3+1,sizeof(unsigned short));
	for (i = 0; i<l1;i++) a1[i]=(i+2)%4+1; a1[i]=0;
	for (i = 0; i<l2;i++) a2[i]=(i+1)%4+1; a2[i]=0;
	for (i = 0; i<l3;i++) a3[i]=(i+1)%4+1; a3[i]=0;
	cout<<"Try"<<endl<<flush;
	try
	{
		sp = new SequencesPile();
    } catch (const DumbException &de) {
        cout << de;
        return 10;
    }
    cout<<"Add"<<endl<<flush;
	try
	{
		sp->add(a1,l1,"seq1");
		sp->add(a2,l2,"seq2");
		sp->add(a3,l3,"seq3");
    } catch (const AtgcException &ae) {
        cout << ae;
        return 10;
    }

    cout<<"Sequence bunch is created..."<<endl<<flush;

	for (unsigned int j=0;j<sp->size();j++)
	{
		cout<<"Lengh of seq "<<j<<"="<<(*sp)[j].size()<<"\n";
		for (unsigned int i=0;i<(*sp)[j].size();cout<<(*sp)[j][i++]<<" ");
		cout<<"\n";
	}
	cout<<sp->size()<<" sequences;  max length is "<<sp->max_length<<endl;

	unsigned short cb=0;
	try
	{
		sc=new SymbolsCounter(*sp,letters,cb);
    } catch (const DumbException &de) {
        cout << de;
        return 10;
    }

    cout<<"Counter was created with pattern length="
						<<pat_len<<" and alphabet="<<letters<<endl<<flush;
	cout<<*sc;
	mcs.positions[0]=0;
	mcs.positions[1]=0;
	mcs.positions[2]=0;
	mcs.SetLength(pat_len, *sp,off);
	cout<<"Calculating.....\n"<<flush;
	sc->calculate(mcs);
	cout<<*sc;
	cout<<"Excluding 1\n"<<flush;
	sc->exclude_sequence(mcs,1);
	cout<<*sc;
	cout<<"Including 1\n"<<flush;
	cout<<*sc;
	sc->include_sequence(mcs,1);
	cout<<"Calculating.....\n"<<flush;
	sc->calculate(mcs);
	cout<<*sc;
	double sum=0;
	long total_let=0;
	cout<<"Backgr:\n";
	for (unsigned int i = 1; i<=letters;i++)
	{
		cout<<sc->background_probability(i)<<"; ";
		sum+=sc->background_probability(i);
	}
	cout<<"\nsum="<<sum<<"\n";
	cout<<"Pseudo:\n";
	sum=0;
	for (unsigned int i = 1; i<=letters;i++)
	{
		cout<<sc->pseudocount(i)<<"; ";
		sum+=sc->pseudocount(i);
	}
	cout<<"\nsum="<<sum<<"\n";
	cout<<"Coded sum="<<sc->pseudocounts_sum()<<"\n";
	cout<<"\nFg:\n";
	for (unsigned int j = 0; j<pat_len;j++)
	{
		sum=0;
		for (unsigned int i = 1; i<=letters;i++)
		{
			cout<<sc->foreground_probability(j,i)<<"; ";
			sum+=sc->foreground_probability(j,i);
		}
		cout<<"\n";
		cout<<"\nsum="<<sum<<"\n";
	}
	cout<<"\n";
	for (unsigned int i = 1; i<=letters;i++)
	{
		cout<<sc->nonpattern_count(i)<<"; ";
		total_let+=sc->nonpattern_count(i);
	}
	cout<<"\nsum of letters="<<total_let<<"\n";
	cout<<"total : "<<sc->all_total()<<"   current total = "<<sc->total();
	cout<<"\n";



	cout<<"Setting pseudoc sum to 10\n";
	sc->change_pseudocounts_sum(10);

	cout<<"Backgr:\n";
	for (unsigned int i = 1; i<=letters;i++)
	{
		cout<<sc->background_probability(i)<<"; ";
		sum+=sc->background_probability(i);
	}
	cout<<"\nsum="<<sum<<"\n";
	cout<<"Pseudo:\n";
	sum=0;
	for (unsigned int i = 1; i<=letters;i++)
	{
		cout<<sc->pseudocount(i)<<"; ";
		sum+=sc->pseudocount(i);
	}
	cout<<"\nsum="<<sum<<"\n";
	cout<<"Coded sum="<<sc->pseudocounts_sum()<<"\n";
	cout<<"\nFg:\n";
	for (unsigned int j = 0; j<pat_len;j++)
	{
		sum=0;
		for (unsigned int i = 1; i<=letters;i++)
		{
			cout<<sc->foreground_probability(j,i)<<"; ";
			sum+=sc->foreground_probability(j,i);
		}
		cout<<"\n";
		cout<<"\nsum="<<sum<<"\n";
	}
	cout<<"\n";
	total_let=0;
	for (unsigned int i = 1; i<=letters;i++)
	{
		cout<<sc->nonpattern_count(i)<<"; ";
		total_let+=sc->nonpattern_count(i);
	}
	cout<<"\nsum of letters="<<total_let<<"\n";
	cout<<"total : "<<sc->all_total()<<"   current total = "<<sc->total();
	cout<<"\n";

	delete sc;
	delete sp;
	cout<<endl;
	return 0;
}

