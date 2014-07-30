[+ AutoGen5 template hpp cpp +][+ CASE (suffix) +][+ ==  hpp  +]
/*  legalistic nonsense:
 *
 * stolen from A. Sheplyakov
[+ (gpl "mssmparam" " *  ") +]
 */

/** \file mssmparam.hpp Declaration of MSSM parameters */
#ifndef _MSSM_MODEL_H
#define _MSSM_MODEL_H
#include <cassert>
#include <stdexcept>
#include <ginac/ginac.h>

namespace MSSM {

/*! \brief helper class to initialize symbols */
class mssm_init {
				static int count;
	public:
				mssm_init();
				~mssm_init();
};

static mssm_init my_initializer;

[+ FOR mssmmass "\n" +]
extern const GiNaC::symbol *_[+ (get "element") +]_p;
extern const GiNaC::symbol &_[+ (get "element") +];
/*! [+ (get "info") +] */
extern const GiNaC::symbol [+ (get "element") +];
[+ ENDFOR mssmmass +]


[+ FOR mssmparam "\n" +]
extern const GiNaC::symbol *_[+ (get "element") +]_p;
extern const GiNaC::symbol &_[+ (get "element") +];
/*! [+ (get "info") +] */
extern const GiNaC::symbol [+ (get "element") +];
[+ ENDFOR mssmparam +]

class io_helper
{
	static int count;
	static GiNaC::lst allsyms;
	static GiNaC::lst allmasses;
		public:
			io_helper();
			inline const GiNaC::lst getall() const { return allsyms; }
			inline const GiNaC::lst get_all_masses() const { return allmasses; }
};

} // namespace MSSM


#endif

[+ ==  cpp  +]
/*  legalistic nonsense:
 *
 *  This file is part of mb2l.
 *
[+ (gpl "mb2l" " *  ") +]
 */

/** \file mssmparam.cpp Initialization of symbols */

#include "mssmparam.hpp"
using namespace std;
using namespace GiNaC;

namespace MSSM {

[+ FOR mssmmass "\n" +]
const symbol *_[+ (get "element") +]_p;
const symbol &_[+ (get "element") +]=*_[+ (get "element") +]_p;
const symbol [+ (get "element") +]=_[+ (get "element") +]; 
[+ ENDFOR mssmmass +]
// total: [+ (count "mssmmass") +] masses


[+ FOR mssmparam "\n" +]
const symbol *_[+ (get "element") +]_p;
const symbol &_[+ (get "element") +]=*_[+ (get "element") +]_p;
const symbol [+ (get "element") +]=_[+ (get "element") +]; 
[+ ENDFOR mssmparam +]

int mssm_init::count=int(0);

mssm_init::mssm_init() {
				if (count++==0) {
	// initialize masses
[+ FOR mssmmass "\n" +]
		_[+ (get "element") +]_p = reinterpret_cast<const symbol*>(&((new symbol("[+ (get "element") +]"))->setflag(status_flags::dynallocated)));
[+ ENDFOR mssmmass +]

	// initialize params
[+ FOR mssmparam "\n" +]
		_[+ (get "element") +]_p = reinterpret_cast<const symbol*>(&((new symbol("[+ (get "element") +]"))->setflag(status_flags::dynallocated)));
[+ ENDFOR mssmparam +]

				}
}

mssm_init::~mssm_init() {
	if (--count==0) {
		//XXX: maybe, need to do some cleanup?
		}
	}

int io_helper::count = int(0);
lst io_helper::allsyms = lst();
lst io_helper::allmasses = lst();

io_helper::io_helper() {
	if (count++==0) {
[+ FOR mssmparam "\n" +]
	allsyms.append([+ (get "element") +]); 
[+ ENDFOR mssmparam +]
[+ FOR mssmmass "\n" +]
	allsyms.append([+ (get "element") +]); 
	allmasses.append([+ (get "element") +]);
[+ ENDFOR mssmmass +]

	}
}

} // namespace MSSM

[+ == m +] [+ COMMENT definitions for Mathematica +]
[+ FOR mssmmass +]
[+ IF (exist? "faName") +]
[+ (get "faName") +] = [+ (get "element") +];[+ ENDIF +]
[+ IF (exist? "real") +]
[+ (get "element") +] /: Conjugate[[+ (get "element")+]] = [+ (get "element") +];
[+ ENDIF +]
m[+ (get "element") +] = Power[[+ (get "element") +], 2];
[+ ENDFOR +]

[+ FOR mssmparam +]
[+ IF (exist? "faName") +]
[+ (get "faName") +] = [+ (get "element") +];
[+ ENDIF +]
[+ IF (or (not (exist? "complex")) (exist? "real")) +]
[+ (get "element") +] /: Conjugate[[+ (get "element")+]] = [+ (get "element") +];
[+ ENDIF +]
[+ ENDFOR +]

If [ $CalcNum,
	[+ FOR mssmmass +]
	[+ IF (exist? "testVal") +][+ (get "element") +] = [+ (get "testVal")+]; [+ ENDIF +]
	[+ ENDFOR +]
	[+ FOR mssmparam +]
	[+ IF (exist? "testVal") +][+ (get "element") +] = [+ (get "testVal")+]; [+ ENDIF +]
	[+ ENDFOR +]

];

[+ ESAC +]

// hey, :vi:ft=cpp



