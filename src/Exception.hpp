/****************************************************************************\
SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2021
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For general describtion of the classes declared in the header, see headers.txt 
$Id$
\****************************************************************************/

#ifndef _EXCEPTION_HPP
#define _EXCEPTION_HPP

#include <iostream>



struct DumbException
{
    std::string info;
    explicit DumbException(const char *str = "")
        : info(str)
    {}
};

inline std::ostream & operator<< (std::ostream & o, const DumbException & de)
{
	o<<de.info;
	return o;
}

struct AtgcException : public DumbException
{
    explicit AtgcException(const char *str = "")
        : DumbException(str)
    {}
};

struct IOStreamException : public DumbException
{
    explicit IOStreamException(const char *str = "")
        : DumbException(str)
    {}
};


#endif //_EXCEPTION_HPP
