#include <iostream>
#include <string>

using namespace std;

#include "../Diagnostics.hpp"


void howtouse()
{
	cout<<"Usage: DiagnosticsTest";
	
};


int main(int argc, char ** argv)
{
	if (argc!=1)
	{
		howtouse();
		return 10;
	}

	Diagnostics diagnostics;

	diagnostics<<"First line.\n";
	diagnostics<<"Second line.\n";
	diagnostics<<"Third line.\n";

    diagnostics.output_mode=txt_output;
	cout<<"Text:\n"<<diagnostics<<endl;
	diagnostics.output_mode=comment_output;
	cout<<"Comment:\n"<<diagnostics<<endl;
	diagnostics.output_mode=html_output;
	cout<<"HTML:\n"<<diagnostics<<endl;
	diagnostics.output_mode=xml_output;
	cout<<"XML:\n"<<diagnostics<<endl;

	return 0;
}
