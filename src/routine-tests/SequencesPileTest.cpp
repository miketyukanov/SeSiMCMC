#include <iostream>
#include <iomanip>

#include "../Sequences.hpp"

int main()
{
	cin.sync_with_stdio();
	SequencesPile sp;
	try{
		cin>>sp;
    } catch (const DumbException &de) {
        cout << de;
    }
    cout << sp.size() << endl;
    for (unsigned int p=0;p<sp.size();p++)
		cout<<sp[p].size()<<" ";
	cout<<endl<<"Min:"<<sp.min_length<<"  Max:"<<sp.max_length<<"  Median:"<<sp.find_median_length()<<endl;
	try {
		cout<<sp;
    } catch (const DumbException &de) {
        cout << de;
    }
    /*	cout<<endl<<"Removing all < 11"<<endl;
	sp.remove_garbage(11);
	cout << sp.size()<<endl;
	for (unsigned int p=0;p<sp.size();p++)
		cout<<sp[p].size()<<" ";
	cout<<endl<<"Min:"<<sp.min_length<<"  Max:"<<sp.max_length<<"  Median:"<<sp.find_median_length()<<endl;
	try {
		cout<<sp;
	} catch(const DumbException &de) {cout<<de;}*/
    return 0;
}
