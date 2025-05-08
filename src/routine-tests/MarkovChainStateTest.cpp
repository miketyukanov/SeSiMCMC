#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

#include "../Exception.hpp"
#include "../MarkovChainState.hpp"
#include "../Sequences.hpp"


void howtouse()
{
	cout<<"Usage: MarkovChainTest --in|--out FastAfilename\n"
			<<"The both params are oblitatory.\n";
	
};


int main(int argc, char ** argv)
{
	//we need 2 parameters (what to do, --in or --out and the FastA);
	SequencesPile sp;
	enum {no,in,out} mode=no;
	string FastAName;
	if (argc==3)
	{
		if (string(argv[1])=="--out")
		{
			mode=out;
			FastAName=argv[2];
		}
		if (string(argv[1])=="--in")
		{
			mode=in;
			FastAName=argv[2];
		}
		if (string(argv[2])=="--out")
		{
			mode=out;
			FastAName=argv[1];
		}
		if (string(argv[2])=="--in")
		{
			mode=in;
			FastAName=argv[1];
		}
	}
	if (mode==no) //parameters error
	{
		howtouse();
		return 10;
	}

	ifstream fastainput(FastAName.c_str());
	if (! fastainput.good()) 
	{
		fprintf(stderr,"Input file %s cannot be opened for read.\n",FastAName.c_str());
		return 10;
	}

	try{
		fastainput>>sp;
    } catch (const DumbException &de) {
        cout << de;
        return 5;
    }

    MarkovChainState state1(sp.size());
	MarkovChainState state2(sp.size());

	if (mode==out)
	{
		state1.SetLength(10,sp,off); state2.SetLength(11,sp,off);
		//state1.SetSymmetricGap(5); state2.SetSymmetricGap(5);
		cout<<state1<<state2;
		return 0;
	}
	if (mode==in)
	{
		try{
		cin>>state1>>state2;
        } catch (const DumbException &de) {
            cerr << de;
            return 1;
        }
        cout<<state1<<state2<<"Scalar product is "<<state2*state1<<endl<<
				"Strict scalar product is "<<StrictScalarProduct(state1,state2)<<endl;
		return 0;
	}
	return 0;
}
