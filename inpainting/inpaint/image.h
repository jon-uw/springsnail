
#ifndef __INC_IMAGE_H__CE89EEA2_2D61_11D3_9988_005004AE2247__
#define __INC_IMAGE_H__CE89EEA2_2D61_11D3_9988_005004AE2247__

#include <windows.h>

//
//	CImage - Fundamental Image Processing Class
//
/////////////////////////////////////////////////////////////////////////////
#define IS_WIN30_DIB(lpbi)  ((*(LPDWORD)(lpbi)) == sizeof(BITMAPINFOHEADER))

class CImage
{

public:

	// Construction / destruction

	CImage();
	CImage(HANDLE hDib);

#ifdef _DEBUG

	// Doesn't return error message, debug only

	CImage(LONG width, LONG height, WORD bitCount);

#endif

	~CImage();

public:

	// Bitmap handle.

	HANDLE	m_hDib;

	// The bitmap is attached, or the bitmap handle should be exported. 
	// CImage destructor doesn't destroyed it. 

	BOOL	m_bIsAttached;

	// Pointer to the BITMAPINFO structure. Variable valid only while memory 
	// locked.

	LPBITMAPINFO	m_pBMI;

	// Pointer to the bitmap bits. Variable valid only while memory locked.

	LPBYTE	m_pPixel;
	

	// The exact bytes for storing one line.

	LONG	m_lWidthBytes;

	LONG    m_BMPwidth;
	LONG    m_BMPheigth;

public:
	void Convert2Gray(CImage* pDestImg);
	BOOL Save(LPCTSTR lpszFilePath,LPBYTE pData);
//	BOOL SaveDCAsFile(HDC hDC, LPCTSTR lpszFilePath);
	BOOL Draw(HDC hDC, POINT ptOriginDst, SIZE sizeDst, RECT rcSrc);
	BOOL Mirror();
	BOOL Negate();
	COLORREF HLStoRGB(const double& H,const double& L,const double& S);
	double HuetoRGB(double m1,double m2,double h);
	void RGBtoHSL(COLORREF rgb,double* H,double* S,double*L);
	BOOL SetPixel(int x,int y,const COLORREF& value);
	int GetBitsPerPixel() const;
	BOOL GetPixel(int x,int y,COLORREF& value);

	// Lock / unlock memory.

	BOOL	Lock();
	BOOL	Unlock();

	// Pointer to bitmap info. Used only when memory is locked.

	inline LPBITMAPINFO	BMI();

	// Pointer to bitmap bits. Used only when memory is locked.

	inline LPBYTE	Pixel();

	// Image dimensions

	inline LONG		Width() const;
	inline LONG		Height() const;

	// Image size in bytes

	inline LONG		Size() const;

	// Actual bytes used in one line.

	inline LONG		WidthBytes() const;

	// 8-bits grayscale pixel manipulations.
	
	inline BYTE		Get(LONG x, LONG y) const;
	inline VOID		Put(LONG x, LONG y, BYTE byPixel);
	inline LPBYTE	Ptr(LONG x, LONG y);

	// 24-bits bitmap pixel I/O.

	inline RGBTRIPLE	GetRGB(LONG x, LONG y) const;
	inline VOID			PutRGB(LONG x, LONG y, RGBTRIPLE rgbtPixel);
	inline RGBTRIPLE*	PtrRGB(LONG x, LONG y);

	// Image create and destory

	BOOL	Create(LONG width, LONG height, WORD bitCount);
	VOID	Release();

	// HANDLE-related functions

	BOOL	Attach(HANDLE hDib);
	HANDLE	Detach();

	// Clear image.

	VOID	Clear(BYTE byInit = 0);

	// Gaussian blur. The source and destinate image can be same. 

	BOOL	Gaussian(CImage* pDest);
	BOOL	GaussianRGB(CImage* pDest);

	// Duplicate image (effective region)

	VOID	Copy(CImage* pDestDIB);

	// Painting routions

	virtual VOID	Draw(HDC hDC, LONG x, LONG y);
	virtual VOID	Stretch(HDC hDC, LONG x, LONG y, LONG width, LONG height);

	// Image I/O.
	
	virtual BOOL	Load(LPCTSTR lpszFilePath);
	virtual BOOL	Save(LPCTSTR lpszFilePath);

	BOOL AddtoFileEnd(FILE *pfile);
	BOOL LoadFromCurPos(FILE *pfile);

	// zdjiang added functions for super-resolution 24-5-2001

	// convert 24 bit bmp to 8 bit gray bmp
	HANDLE ConvertToGray();

	// extract r, g, and b channel
	HANDLE ExtractChannel(char cFlag);

	// Zoom in using bilinear interplolation
	HANDLE BilinearZoomIn(int ratio);
	double GetPSNRColor(CImage* pRecon,int iBW,int iBH); // 24bits bmp
	double GetPSNRGray(CImage* pRecon,int iBW,int iBH); // 24bits bmp
	CImage* TransToYIQ();

private:
	CImage* CombineYIQ(float* pYData,float* pIData, float* pQData);
	BOOL GetYIQ(float* pYData,float* pIData, float* pQData);
	

protected:
	LONG	GetColorEntries(LPBITMAPINFO pBMI);
};

inline LONG CImage::Width() const
{
	return m_pBMI->bmiHeader.biWidth;
}

/////////////////////////////////////////////////////////////////////////////

inline LPBITMAPINFO CImage::BMI()
{
	return m_pBMI;
}

inline LPBYTE CImage::Pixel()
{
	return m_pPixel;
}

/////////////////////////////////////////////////////////////////////////////

inline LONG CImage::Height() const
{
	return m_pBMI->bmiHeader.biHeight;
}

inline LONG CImage::Size() const
{
	return m_pBMI->bmiHeader.biSizeImage;
}

inline LONG CImage::WidthBytes() const
{
	return m_lWidthBytes;
}

/////////////////////////////////////////////////////////////////////////////

inline BYTE CImage::Get(LONG x, LONG y) const
{
	return *(m_pPixel + y * m_lWidthBytes + x);
}

inline void CImage::Put(LONG x, LONG y, BYTE byPixel)
{
	*(m_pPixel + y * m_lWidthBytes + x) = byPixel;
}

inline LPBYTE CImage::Ptr(LONG x, LONG y)
{
	return m_pPixel + y * m_lWidthBytes + x;
}

/////////////////////////////////////////////////////////////////////////////

inline RGBTRIPLE CImage::GetRGB(LONG x, LONG y) const
{
	return *(RGBTRIPLE*) (m_pPixel + y * m_lWidthBytes + 3 * x);
}

inline void CImage::PutRGB(LONG x, LONG y, RGBTRIPLE rgbtPixel)
{
	*((RGBTRIPLE*)(m_pPixel + y * m_lWidthBytes + 3 * x)) = rgbtPixel;
}

inline RGBTRIPLE* CImage::PtrRGB(LONG x, LONG y)
{
	return (RGBTRIPLE*) (m_pPixel + y * m_lWidthBytes + 3 * x);
}


#endif	// __INC_IMAGE_H__CE89EEA2_2D61_11D3_9988_005004AE2247__
