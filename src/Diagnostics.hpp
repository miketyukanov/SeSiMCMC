/****************************************************************************\
SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2007
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For general describtion of the classes declared in the header, see headers.txt
$Id$
\****************************************************************************/

#ifndef _DIAGNOSTICS_HPP_
#define _DIAGNOSTICS_HPP_

#include <sstream>

enum output_mode_type {
    unknown_output_mode = 0,
    txt_output,
    comment_output,
    html_output,
    xml_output
};

enum status_type { OK, unreliable, fatal };

class Diagnostics : public std::ostringstream
{
	Diagnostics(const Diagnostics &);
	Diagnostics & operator=(const Diagnostics &);
public:
    Diagnostics()
        : std::ostringstream(std::ostringstream::out)
        , output_mode(txt_output)
        , status(OK)
    {}
    output_mode_type output_mode;
    status_type status;
    void do_text_output(std::ostream &o) const;
    void do_comment_output(std::ostream &o) const;
    void do_html_output(std::ostream &o) const;
    void do_xml_output(std::ostream &o) const;
};

inline std::ostream &operator<<(std::ostream &o, const Diagnostics &d)
{
	switch (d.output_mode)
	{
	case txt_output: d.do_text_output(o);break;
	case comment_output: d.do_comment_output(o);break;
	case html_output:d.do_html_output(o);break;
	case xml_output:d.do_xml_output(o);break;
	default: d.do_text_output(o);break;
	}
	return o;
}

inline void Diagnostics::do_text_output(std::ostream &o) const
{
    auto str = std::ostringstream::str();
    for (const auto &val : str)
        o << val;
}

inline void Diagnostics::do_comment_output(std::ostream &o) const
{
	o<<"#";
    auto str = std::ostringstream::str();
    for(const auto& val: str)
	{
		o<<val;
		if (val == '\n') o<<"#";
	}
}

inline void Diagnostics::do_html_output(std::ostream &o) const
{
	o<<"<TT>"<<endl;
    auto str = std::ostringstream::str();
    for(const auto& val: str)
	{
		o<<val;
		if (val == '\n') o<<"<br>\n";
	}
	o<<"</TT>"<<endl;
}

inline void Diagnostics::do_xml_output(std::ostream &o) const
{
	o<<"<comment name=\"diagnostics\" ";
	switch (status)
	{
	case OK: o<<" status=\"OK\"";break;
	case unreliable: o<<" status=\"unreliable\"";break;
	case fatal: o<<" status=\"fatal\"";break;
	}
	o<<">"<<endl;
    auto str = std::ostringstream::str();
    for(const auto& val: str) o<<val;
	o<<"</comment>"<<endl;
}
#endif /*_DIAGNOSTICS_HPP_*/
