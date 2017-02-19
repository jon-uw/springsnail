/*Author: Qiushuang Zhang */
/*E-mail: qszhang@cc.gatech.edu */
/*Nov.29, 2005 */


#include "StdAfx.h"
#include ".\inpainting.h"
#include <memory.h>
#include <math.h>
#define MAX(a, b)  (((a) > (b)) ? (a) : (b)) 
#define MIN(a, b)  (((a) < (b)) ? (a) : (b)) 

inpainting::inpainting(char * name)
{
	COLORREF temp;
	m_pImage = new CImage;
	m_bOpen = m_pImage->Load(name);  //load the image
	if(m_bOpen)
		{	
			m_pImage->Lock();		
			m_width = m_pImage->m_BMPwidth;
			m_height = m_pImage->m_BMPheigth;

			m_mark = new int[m_width*m_height];
			m_confid = new double[m_width*m_height];
			memset( m_confid, 0, m_width*m_height*sizeof(double) );
			m_pri = new double[m_width*m_height];
			m_gray  = new double[m_width*m_height];
			m_source = new bool[m_width*m_height];
			m_color = new COLORREF[m_width*m_height];
			m_r = new double[m_width*m_height];
			m_g = new double[m_width*m_height];
			m_b = new double[m_width*m_height];
		}
	else printf("one or more file is not opened!\n");
	m_top = m_height;  // initialize the rectangle area
    m_bottom = 0; 
	m_left = m_width;
	m_right = 0;
	
	for(int y = 0; y<m_height; y++)
		for(int x = 0; x<m_width; x++)
		{
			m_pImage->GetPixel(x,y,temp);
			m_color[y*m_width+x] = temp;
			m_r[y*m_width+x] = GetRValue(temp);
			m_g[y*m_width+x] = GetGValue(temp);
			m_b[y*m_width+x] = GetBValue(temp);
		}

}

inpainting::~inpainting(void)
{
	if(m_pImage)delete m_pImage;
	if(m_mark)delete m_mark;
	if(m_source)delete m_source;
	if(m_color)delete m_color;
	if(m_r)delete m_r;
	if(m_g)delete m_g;
	if(m_b)delete m_b;
}


void inpainting::Convert2Gray(void)
{
	COLORREF  cc;
	double r,g,b;
	for(int y = 0; y<m_height; y++)
		for(int x = 0; x<m_width; x++)
		{
			m_pImage->GetPixel(x,y,cc);
			r = GetRValue(cc);
			g = GetGValue(cc); 
			b = GetBValue(cc);
			m_gray[y*m_width+x] = (double)((r*3735 + g*19267 + b*9765)/32767);
		}
}

bool inpainting::process(void)
{
	char path[200];
	char temp[30];
	Convert2Gray();  // convert it to gray image
	DrawBoundary();  // first time draw boundary
	draw_source();   // find the patches that can be used as sample texture
	memset(m_pri, 0, m_width*m_height*sizeof(double));
	for(int j= m_top; j<=m_bottom; j++)
	    for(int i = m_left; i<= m_right; i++)if(m_mark[j*m_width+i] == BOUNDARY)m_pri[j*m_width+i] = priority(i,j);//if it is boundary, calculate the priority
	int count=0;
	while(TargetExist()) 
	{
		count++;
		double max_pri = 0;
		int pri_x,pri_y;
		for(int j= m_top; j<=m_bottom; j++)
	    for(int i = m_left; i<= m_right; i++)
		{
			if(m_mark[j*m_width+i] == BOUNDARY&&m_pri[j*m_width+i]>max_pri)// find the boundary pixel with highest priority
			{
				pri_x = i;
				pri_y = j;
				max_pri = m_pri[j*m_width+i];
			}
		}
		//printf("pri_x is %d, pri_y is %d, amount is %f\n", pri_x, pri_y, max_pri);
		int patch_x, patch_y;
		PatchTexture(pri_x, pri_y, patch_x, patch_y);  // find the most similar source patch
		update(pri_x, pri_y, patch_x,patch_y, ComputeConfidence(pri_x,pri_y));// inpaint this area and update confidence
		UpdateBoundary(pri_x, pri_y); // update boundary near the changed area
		UpdatePri(pri_x, pri_y);  //  update priority near the changed area
	
		//m_pImage->Save("inpaint_result.bmp");
		strcpy(path, save_path);
		if(count<10)strcat(path,"000");
		else if(count<100)strcat(path,"00");
		else if(count<1000)strcat(path,"0");
		_itoa( count, temp, 10 );
		strcat(path,temp);
		strcat(path, ".bmp");
		m_pImage->Save(path);
	}
	return true;
}

void inpainting::DrawBoundary(void)// fist time draw boundary
{
	COLORREF color;	
	for(int y = 0; y<m_height; y++)
		for(int x = 0; x<m_width; x++)
		{
			m_pImage->GetPixel(x,y,color);   
			if(color==PAINT_COLOR)// if the pixel is specified as boundary
			{
				m_mark[y*m_width+x] = TARGET;
				m_confid[y*m_width+x] = 0;
			}
			else {
				m_mark[y*m_width+x] = SOURCE;
				m_confid[y*m_width+x] = 1;
			}
		}

	for(int j= 0; j< m_height; j++)
	    for(int i = 0; i< m_width; i++)
		{
			if(m_mark[j*m_width+i]==TARGET)
			{
				if(i<m_left)m_left = i; // resize the rectangle to the range of target area
				if(i>m_right)m_right = i;
				if(j>m_bottom)m_bottom = j;
				if(j<m_top)m_top = j;
				//if one of the four neighbours is source pixel, then this should be a boundary
				if(j==m_height-1||j==0||i==0||i==m_width-1||m_mark[(j-1)*m_width+i]==SOURCE||m_mark[j*m_width+i-1]==SOURCE
					||m_mark[j*m_width+i+1]==SOURCE||m_mark[(j+1)*m_width+i]==SOURCE)m_mark[j*m_width+i] = BOUNDARY;
			}
		}
}



double inpainting::priority(int i, int j)
{
	double confidence, data;
	confidence = ComputeConfidence(i,j); // confidence term
	data = ComputeData(i,j);// data term
	return confidence*data;
}

double inpainting::ComputeConfidence(int i, int j) 
{
	double confidence=0;
	for(int y = MAX(j -winsize,0); y<= MIN(j+winsize,m_height-1); y++)
		for(int x = MAX(i-winsize,0); x<=MIN(i+winsize, m_width-1); x++)
			confidence+= m_confid[y*m_width+x];
	confidence/= (winsize*2+1)*(winsize*2+1);
	return confidence;

}
double inpainting::ComputeData(int i, int j)
{
	gradient grad, temp, grad_T;
	grad.grad_x=0;
	grad.grad_y=0;
	double result;
	double magnitude;
	double max=0;
	int x, y;
	for(y = MAX(j -winsize,0); y<= MIN(j+winsize,m_height-1); y++)
		for( x = MAX(i-winsize,0); x<=MIN(i+winsize, m_width-1); x++)
		{
			// find the greatest gradient in this patch, this will be the gradient of this pixel(according to "detail paper")
			if(m_mark[y*m_width+x] >=0) // source pixel
			{
				//since I use four neighbors to calculate the gradient, make sure this four neighbors do not touch target region(big jump in gradient)
				if(m_mark[y*m_width+x+1]<0||m_mark[y*m_width+x-1]<0||m_mark[(y+1)*m_width+x]<0||m_mark[(y-1)*m_width+x]<0)continue;
 				temp = GetGradient(x,y); 
				magnitude = temp.grad_x*temp.grad_x+temp.grad_y*temp.grad_y;
				if(magnitude>max)
				{
					grad.grad_x = temp.grad_x;
					grad.grad_y = temp.grad_y;
					max = magnitude;
				}
			}
		}
		grad_T.grad_x = grad.grad_y;// perpendicular to the gradient: (x,y)->(y, -x)
		grad_T.grad_y = -grad.grad_x;

	norm nn = GetNorm(i,j);
	result = nn.norm_x*grad_T.grad_x+nn.norm_y*grad_T.grad_y; // dot product
	result/=255; //"alpha" in the paper: normalization factor
	result = fabs(result);			
	return result;
}



gradient inpainting::GetGradient(int i, int j)
{
	gradient result;
	result.grad_x = (m_gray[j*m_width+i+1] - m_gray[j*m_width+i-1])/2.0;
	result.grad_y = (m_gray[(j+1)*m_width +i] - m_gray[(j-1)*m_width+i])/2.0;

	if(i==0)result.grad_x = m_gray[j*m_width+i+1] - m_gray[j*m_width+i];
	if(i==m_width-1)result.grad_x = m_gray[j*m_width+i] - m_gray[j*m_width+i-1];
	if(j==0)result.grad_y = m_gray[(j+1)*m_width +i] - m_gray[j*m_width+i];
	if(j==m_height-1)result.grad_y = m_gray[j*m_width +i] - m_gray[(j-1)*m_width+i];
	return result;
}

norm inpainting::GetNorm(int i, int j)
{
	norm result;
	int num=0;
	int neighbor_x[9];
	int neighbor_y[9];
	int record[9];
	int count = 0;
	for(int y = MAX(j-1,0); y<=MIN(j+1,m_height-1); y++)
		for(int x = MAX(i-1,0); x<=MIN(i+1,m_width-1); x++)
		{
			count++;
			if(x==i&&y==j)continue;
			if(m_mark[y*m_width+x]==BOUNDARY)
			{
				num++;
				neighbor_x[num] = x;
				neighbor_y[num] = y;		
				record[num]=count;
			}
		}
		if(num==0||num==1) // if it doesn't have two neighbors, give it a random number to proceed
		{
			result.norm_x = 0.6;
			result.norm_y = 0.8;
			return result;
		}
		// draw a line between the two neighbors of the boundary pixel, then the norm is the perpendicular to the line
			int n_x = neighbor_x[2]-neighbor_x[1];
			int n_y = neighbor_y[2]-neighbor_y[1];
			int temp=n_x;
			n_x = n_y;
			n_y = temp;
			double square = pow(double(n_x*n_x + n_y*n_y),0.5);
		
	result.norm_x = n_x/square;
	result.norm_y =n_y/square;
	return result;
}

bool inpainting::draw_source(void) 
{
	// draw a window around the pixel, if all of the points within the window are source pixels, then this patch can be used as a source patch
	bool flag;
		for(int j = 0; j<m_height; j++)
		for(int i = 0; i<m_width; i++)
		{
			flag=1;
			if(i<winsize||j<winsize||i>=m_width-winsize||j>=m_height-winsize)m_source[j*m_width+i]=false;//cannot form a complete window
			else
			{
				for(int y = j-winsize; y<=j+winsize; y++)
				{
				for(int x = i-winsize; x<=i+winsize; x++)
				{
					if(m_mark[y*m_width+x]!=SOURCE)
					{
						m_source[j*m_width+i]=false;
						flag = false;
						break;			
					}		
				}
				if(flag==false)break;
				}
			    if(flag!=false)m_source[j*m_width+i]=true;
			}
		}
	return true;
}

bool inpainting::PatchTexture(int x, int y, int &patch_x, int &patch_y)
{
	double temp_r;
	double temp_g;
	double temp_b;
	// find the most similar patch, according to SSD
		COLORREF color_target, color_source, color_diff;
		double r0,g0,b0;
		double r1,g1,b1;
	
		long min=99999999;
		long sum;
		int source_x, source_y;
		int target_x, target_y;
		for(int j = 0; j<m_height; j++)
		for(int i = 0; i<m_width; i++)
		{		
			if(m_source[j*m_width+i]==false)continue;
			sum=0;
			for(int iter_y=(-1)*winsize; iter_y<=winsize; iter_y++)
				for(int iter_x=(-1)*winsize; iter_x<=winsize; iter_x++)
				{
					source_x = i+iter_x;
					source_y = j+iter_y;

					target_x = x+iter_x;
					target_y = y+iter_y;
					
					if(target_x<0||target_x>=m_width||target_y<0||target_y>=m_height)continue;

					if(m_mark[target_y*m_width+target_x]>=0)
					{
						temp_r = m_r[target_y*m_width+target_x]-m_r[source_y*m_width+source_x];
						temp_g = m_g[target_y*m_width+target_x]-m_g[source_y*m_width+source_x];
						temp_b = m_b[target_y*m_width+target_x]-m_b[source_y*m_width+source_x];
						sum+= temp_r*temp_r + temp_g*temp_g + temp_b*temp_b;
					}
				}
				if(sum<min)
				{
					min=sum;
					patch_x = i;
					patch_y = j;
				}
		}
	//	printf("patch_x is %d, patch_y is %d\n", patch_x, patch_y);
	return true;
}

bool inpainting::update(int target_x, int target_y, int source_x, int source_y, double confid)
{
	COLORREF color;
	double r,g,b;

	int x0,y0,x1,y1;
	for(int iter_y=(-1)*winsize; iter_y<=winsize; iter_y++)
		for(int iter_x=(-1)*winsize; iter_x<=winsize; iter_x++)
		{
			x0 = source_x+iter_x;
			y0 = source_y+iter_y;

			x1 = target_x + iter_x;
			y1 = target_y + iter_y;
			if(m_mark[y1*m_width+x1]<0)
			{
				m_pImage->SetPixel(x1,y1,m_color[y0*m_width+x0]);// inpaint the color
				m_color[y1*m_width+x1] = m_color[y0*m_width+x0];
				m_r[y1*m_width+x1] = m_r[y0*m_width+x0];
				m_g[y1*m_width+x1] = m_g[y0*m_width+x0];
				m_b[y1*m_width+x1] = m_b[y0*m_width+x0];
			m_gray[y1*m_width+x1] = (double)((m_r[y0*m_width+x0]*3735 + m_g[y0*m_width+x0]*19267 + m_b[y0*m_width+x0]*9765)/32767); // update gray image
				m_confid[y1*m_width+x1] = confid; // update the confidence
			}
		}
	return true;
}

bool inpainting::TargetExist(void)
{
		for(int j= m_top; j<=m_bottom; j++)
	    for(int i = m_left; i<= m_right; i++)if(m_mark[j*m_width+i]<0)return true;
	return false;
}

void inpainting::UpdateBoundary(int i, int j)// just update the area near the changed patch. (+-2 pixels)
{
	int x, y;
	COLORREF color;	

	for(y = MAX(j -winsize-2,0); y<= MIN(j+winsize+2,m_height-1); y++)
		for( x = MAX(i-winsize-2,0); x<=MIN(i+winsize+2, m_width-1); x++)
		{
			m_pImage->GetPixel(x,y,color);   // if the pixel is specified as boundary
			if(color==PAINT_COLOR)m_mark[y*m_width+x] = TARGET;
			else m_mark[y*m_width+x] = SOURCE;
		}

	for(y = MAX(j -winsize-2,0); y<= MIN(j+winsize+2,m_height-1); y++)
		for( x = MAX(i-winsize-2,0); x<=MIN(i+winsize+2, m_width-1); x++)
		{
			if(m_mark[y*m_width+x]==TARGET)
			{
				if(y==m_height-1||y==0||x==0||x==m_width-1||m_mark[(y-1)*m_width+x]==SOURCE||m_mark[y*m_width+x-1]==SOURCE
					||m_mark[y*m_width+x+1]==SOURCE||m_mark[(y+1)*m_width+x]==SOURCE)m_mark[y*m_width+x] = BOUNDARY;
			}
		}		
}

void inpainting::UpdatePri(int i, int j) // just update the area near the changed patch. (+-3 pixels)
{
	int x,y;
	for(y = MAX(j -winsize-3,0); y<= MIN(j+winsize+3,m_height-1); y++)
		for( x = MAX(i-winsize-3,0); x<=MIN(i+winsize+3, m_width-1); x++)if(m_mark[y*m_width+x] == BOUNDARY)m_pri[y*m_width+x] = priority(x,y);

}
