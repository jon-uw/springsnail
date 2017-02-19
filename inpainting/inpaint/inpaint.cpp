// inpaint.cpp : Defines the entry point for the console application.
//
/*Author: Qiushuang Zhang */
/*E-mail: qszhang@cc.gatech.edu */
/*Nov.29, 2005 */

#include "stdafx.h"
#include "image.h"
#include "inpainting.h"


int _tmain(int argc, _TCHAR* argv[])
{
	inpainting test("test.bmp");
	test.process();
	
	return 0;
}

