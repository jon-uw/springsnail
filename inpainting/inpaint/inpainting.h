#pragma once

/*Author: Qiushuang Zhang */
/*E-mail: qszhang@cc.gatech.edu */
/*Nov.29, 2005 */

#ifndef INPAINTING_H
#define INPAINTING_H

#include "image.h"
#define save_path ".\\result"

#define PAINT_COLOR  RGB(0,255,0)  // the color used to specify the target area
#define SOURCE 0
#define TARGET -1
#define BOUNDARY -2
#define winsize 4  // the window size

typedef struct
{
	double grad_x;
	double grad_y;
}gradient; //the structure that record the gradient

typedef struct
{
	double norm_x;
	double norm_y;
}norm;  // the structure that record the norm

class inpainting
{
public:
	
	CImage * m_pImage; // the image to be inpainted
	bool m_bOpen; // whether it is successfully opened
	int m_width; // image width
	int m_height; // image height

	COLORREF * m_color;
	double * m_r;
	double * m_g;
	double * m_b;
	
	int m_top, m_bottom, m_left, m_right; // the rectangle of inpaint area

	int * m_mark;// mark it as source(0) or to-be-inpainted area(-1) or bondary(-2).
	double * m_confid;// record the confidence for every pixel
	double * m_pri; // record the priority for pixels. only boudary pixels will be used
	double * m_gray; // the gray image
	bool * m_source; // whether this pixel can be used as an example texture center

	inpainting(char * name);
	~inpainting(void);
	bool process(void);  // the main function to process the whole image
	void DrawBoundary(void);  // the first time to draw boundary on the image.

	double ComputeConfidence(int i, int j); // the function to compute confidence
	double priority(int x, int y); // the function to compute priority
	double ComputeData(int i, int j);//the function to compute data item
	void Convert2Gray(void);  // convert the input image to gray image
	gradient GetGradient(int i, int j); // calculate the gradient at one pixel
	norm GetNorm(int i, int j);  // calculate the norm at one pixel
	bool draw_source(void);  // find out all the pixels that can be used as an example texture center
	bool PatchTexture(int x, int y,int &patch_x,int &patch_y);// find the most similar patch from sources.
	bool update(int target_x, int target_y, int source_x, int source_y, double confid);// inpaint this patch and update pixels' confidence within this area
	bool TargetExist(void);// test whether this is still some area to be inpainted.
	void UpdateBoundary(int i, int j);// update boundary
	void UpdatePri(int i, int j); //update priority for boundary pixels.
};


#endif
