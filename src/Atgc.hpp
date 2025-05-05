/*****************************************************************************\
SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2021
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For general describtion of the classes declared in the header, see headers.txt 
$Id$
\*****************************************************************************/

#ifndef _ATGC_HPP
#define _ATGC_HPP

//#include <ctype.h>

#include <iostream>
#include <iomanip>

#include <vector>
#include <string>

#include <stdio.h>
#include <math.h>


#include "Exception.hpp"

extern "C"
{
	#include "Random.h"
}
//Here, we define unsigned short <-> atgc letter functions and
//procedure to output unsigned short * array faragment as 
//atgc sequence fragment.
//
//We do not define a definite class line 
//AtgcVector::public vector<unsigned short> for we are not sure what it need

const unsigned short MaxAtgcSymbol = 4;

const char atgc[5]={'x','a','t','g','c'};

//s = g || c (strong)
//w = a || t (weak)

//r = a || g (purines)
//y = c || t (pyrimidines)

//k = g || t
//m = a || c

//b = g || t || c (no - a)
//h = a || t || c (no - g)

//n = any residue, now the only used.
// a set of 1,2,3,4 correspods to a,t,g,c 

class Atgc
{
public:
	static char ushort2atgc(unsigned short symbol);
	static unsigned short atgc2ushort(char letter);
	static unsigned short complement(unsigned short sym);
	static char complement(char sym);

    static std::vector<unsigned short> &complement(const std::vector<unsigned short> &source,
                                                   std::vector<unsigned short> &dest);

    static std::vector<unsigned short> &complement(
        std::vector<unsigned short> &source); //in-place version

    static std::string &complement(const std::string &source, std::string &dest);

    static std::string &complement(std::string &dest); //in-place version

    static std::string &atgc2string(const std::vector<unsigned short> &source,
                                    std::string &dest = *new(std::string));

    static std::vector<unsigned short> &string2atgc(
        const std::string &source,
        std::vector<unsigned short> &dest = *new(std::vector<unsigned short>));
};

inline
unsigned short Atgc::atgc2ushort(char letter) 
{
	char a=tolower(letter);
	if (a=='a') return 1;
	if (a=='t') return 2;
	if (a=='g') return 3;
	if (a=='c') return 4;
	//we suppose that the random generator is already inited.
	if (a=='n') 
	{
		int r=(int)(floorf(4.*uni()))+1;
		if (r==5) r=4;
		return r;
	};
		
	if (a=='s')
	{
		if (uni()>.5) return 3;
		return 4;
	}
	
	if (a=='w')
	{
		if (uni()>.5) return 1;
		return 2;
	}
	
	if (a=='r')
	{
		if (uni()>.5) return 1;
		return 3;
	}
	if (a=='y')
	{
		if (uni()>.5) return 2;
		return 4;
	}
	if (a=='x') //it is for the additional mask symbol 'x' 
	{
		//int r=(int)(floorf(4.*uni()))+1;
		//if (r==5) r=4;
		return 0;
		//we cannot define the masked symbol,
		//whatever.
	};
	if (a=='-') //it is for the additional mask symbol '-' 
	{
		//int r=(int)(floorf(4.*uni()))+1;
		//if (r==5) r=4;
		return 0;
		//we cannot define the masked symbol,
		//whatever.
	};
    std::string message = "Trying to read symbol \'";
    message+=letter;
	message+="\' as nucleoutide.\n";
	throw * new AtgcException(message.c_str());
	return 0;
}

inline
char Atgc::ushort2atgc(unsigned short symbol)
{
	if (symbol>MaxAtgcSymbol)
	{
        std::string message = "Trying to interpret ";
        char symb_no[10];
	  snprintf(symb_no,9,"%1i", symbol);	
		message+=symb_no;
		message+=" as an atgc number (1..4).\n";
		throw * new AtgcException("Trying to get atgc char from something other then 1..4.\n");
	}
	return atgc[symbol];
}

inline
unsigned short Atgc::complement(unsigned short sym)
{
//	1<->2
//	3<->4
//	0<->0
		return sym?sym+(sym%2)*2-1:0;
}

inline
char Atgc::complement(char sym)
{
	unsigned short up=isupper(sym);
	char a=tolower(sym);
	if (a=='a') return up?toupper('t'):'t';
	if (a=='t') return up?toupper('a'):'a';
	if (a=='g') return up?toupper('c'):'c';
	if (a=='c') return up?toupper('g'):'g';
	if (a=='n') return up?toupper('n'):'n'; 
	if (a=='w') return up?toupper('w'):'w';
	if (a=='s') return up?toupper('s'):'s';
	if (a=='r') return up?toupper('y'):'y';
	if (a=='y') return up?toupper('r'):'r';
	if (a=='x') return up?toupper('x'):'x';
    std::string message = "Trying to fing complement to symbol \'";
    message+=sym;
	message+="\'.\n";
	throw * new AtgcException(message.c_str());
	return 0;
}

inline std::vector<unsigned short> &Atgc::complement(std::vector<unsigned short> &dest)
{
	//we do not use STL swap because it is more quick to swap and 
	//reverse bases in one pass.
	//on the other hand, we want the algoryth to be stable for in-place
	//operations. the part is not time_critical.
    std::vector<unsigned short>::iterator l = dest.begin();
    std::vector<unsigned short>::iterator r = dest.end();
    r--;
    while (l<=r)
	{
		unsigned short buf=*l;
		*l++=complement(*r);
		*r--=complement(buf);
	}
	return dest;
}

inline std::vector<unsigned short> &Atgc::complement(const std::vector<unsigned short> &source,
                                                     std::vector<unsigned short> &dest)
{
	//we do not use STL swap because it is more quick to swap and 
	//reverse bases in one pass.
	//on the other hand, we want the algoryth to be stable for in-place
	//operations. the part is not time_critical.
	if (&source==&dest) return complement(dest);
	dest.clear();
    std::vector<unsigned short>::const_iterator r = source.end();
    r--;
    while (source.begin()<=r)
		dest.push_back(complement(*r--));
	return dest;
}

inline std::string &Atgc::complement(std::string &dest) //in-place version
{
	//we do not use STL swap because it is more quick to swap and 
	//reverse bases in one pass.
	//on the other hand, we want the algoryth to be stable for in-place
	//operations. the part is not time_critical.
    std::string::iterator l = dest.begin();
    std::string::iterator r = dest.end();
    r--;
	while (l<=r)
	{
		char buf=*l;
		*l++=complement(*r);
		*r--=complement(buf);
	}
	return dest;
}

inline std::string &Atgc::complement(const std::string &source, std::string &dest)
{
	//we do not use STL swap because it is more quick to swap and 
	//reverse bases in one pass.
	//on the other hand, we want the algoryth to be stable for in-place
	//operations. the part is not time_critical.
	if (&source==&dest) return complement(dest);
	dest.erase();
    std::string::const_iterator r = source.end();
    r--;
	while (source.begin()<=r)
		dest.push_back(complement(*r--));
	return dest;
}

inline std::string &Atgc::atgc2string(const std::vector<unsigned short> &source, std::string &dest)
{
    dest="";
    std::vector<unsigned short>::const_iterator r = source.begin();
    while (r != source.end())
        dest.push_back(ushort2atgc(*r++));
	return dest;
}

inline std::vector<unsigned short> &Atgc::string2atgc(const std::string &source,
                                                      std::vector<unsigned short> &dest)
{
    dest.clear();
    std::string::const_iterator r = source.begin();
    while ( r!=source.end() )
		dest.push_back(atgc2ushort(*r++));
	return dest;
}

#endif // _ATGC_HPP
