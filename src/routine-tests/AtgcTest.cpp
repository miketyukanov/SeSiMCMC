#include <vector>
#include <string>

#include "../Atgc.hpp"

int main()
{
    std::vector<unsigned short> atgc,atgc1;
    std::string str;
	str="attttaaaagcccc";
	Atgc::string2atgc(str,atgc);
    std::cout<<"Original string:\n"<<str<<std::endl
            <<"Back converted:"<<Atgc::atgc2string(atgc,str)<<std::endl;
	Atgc::complement(atgc,atgc1);
    std::cout<<"Copy complement:"<<Atgc::atgc2string(atgc1,str)<<std::endl;
	Atgc::complement(atgc);
    std::cout<<"Inplace complement:"<<Atgc::atgc2string(atgc,str)<<std::endl;
	Atgc::complement(atgc,atgc);
    std::cout<<"One more inplace :"<<Atgc::atgc2string(atgc,str)<<std::endl;
}
