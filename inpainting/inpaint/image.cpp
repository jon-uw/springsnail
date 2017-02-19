 /////////////////////////////////////////////////////////////////////////////
//
//	Image . cpp
//
//	8 / 24 bits Device Independent Bitmap manipulations
//
//	Copyright (c) 1999 by ArcSoft Inc.
//
//	Version 2.0
//
//	Last modified: June 1999
//
//	Programmer: Alex Zhou
//
/////////////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "image.h"


/////////////////////////////////////////////////////////////////////////////
//
//	class CImage
//
/////////////////////////////////////////////////////////////////////////////

CImage::CImage()
{
	m_bIsAttached = FALSE;

	m_hDib = NULL;

	m_pBMI = NULL;
	m_pPixel = NULL;
}

CImage::CImage(HANDLE hDib)
{
	CImage();
	m_bIsAttached = TRUE;

	Attach(hDib);
}

#ifdef _DEBUG

CImage::CImage(LONG width, LONG height, WORD bitCount)
{
	CImage();

	if (!Create(width, height, bitCount) )
		MessageBox(NULL, "Memory Insufficient.", "CImage Run-Time Error Message", MB_OK);
}

#endif

CImage::~CImage()
{
	Release();
}

/////////////////////////////////////////////////////////////////////////////

BOOL CImage::Lock()
{
	if ( ! ( m_pBMI = (LPBITMAPINFO) GlobalLock(m_hDib) ) )
		return FALSE;

	switch (m_pBMI->bmiHeader.biBitCount)
	{
	case 1:
		m_pPixel = (LPBYTE) m_pBMI + sizeof(BITMAPINFOHEADER) + 2 * sizeof(RGBQUAD);
		break;

	case 8:
		m_pPixel = (LPBYTE) m_pBMI + sizeof(BITMAPINFOHEADER) + 256 * sizeof(RGBQUAD);
		break;
	case 16:
	case 24:
		m_pPixel = (LPBYTE) m_pBMI + sizeof(BITMAPINFOHEADER);
		break;

	default:
		assert(FALSE);
	}
	

	return TRUE;
}

BOOL CImage::Unlock()
{
	return GlobalUnlock(m_hDib);
}

/////////////////////////////////////////////////////////////////////////////
	
LONG CImage::GetColorEntries(LPBITMAPINFO pBMI)
{
	assert(pBMI);

	if (pBMI->bmiHeader.biClrUsed != 0)
		return pBMI->bmiHeader.biClrUsed;

	switch (pBMI->bmiHeader.biBitCount)
	{
	
	case 1:
		return 2;

	case 4:
		return 16;
	
	case 8:
		return 256;
	
	default:
		return 0;
	
	}
}

// zdjiang added functions for super-resolution 24-5-2001

HANDLE CImage::BilinearZoomIn(int ratio)
{
	assert(m_pBMI->bmiHeader.biBitCount == 8 || m_pBMI->bmiHeader.biBitCount == 24);

	int nZoomWidth = Width() * ratio;
	int nZoomHeight = Height() * ratio;
	int bitCount = m_pBMI->bmiHeader.biBitCount;
	
	HANDLE hDib;
	LPBITMAPINFO pBmpInfo;
	LONG lColors, lColorTableSize, lBitsSize, lBISize, lWidthBytes;

	BITMAPINFOHEADER BmpInfoHdr;

	BmpInfoHdr.biSize = sizeof(BITMAPINFOHEADER);
	BmpInfoHdr.biWidth = nZoomWidth;
	BmpInfoHdr.biHeight = nZoomHeight;
	BmpInfoHdr.biPlanes = 1;
	BmpInfoHdr.biBitCount = bitCount;
	BmpInfoHdr.biCompression = BI_RGB;
	BmpInfoHdr.biSizeImage = (nZoomWidth * bitCount + 31) / 32 * 4 * nZoomHeight;
	BmpInfoHdr.biXPelsPerMeter = 0;
	BmpInfoHdr.biYPelsPerMeter = 0;
	BmpInfoHdr.biClrUsed = 0;
	BmpInfoHdr.biClrImportant = 0;

	// The number of bytes for one scan line.

	lWidthBytes = (BmpInfoHdr.biWidth * BmpInfoHdr.biBitCount + 31) / 32 * 4;

	// Determine memory needed.
	switch (BmpInfoHdr.biBitCount)
	{
	case 1:
		lColors = 2;
		BmpInfoHdr.biClrUsed = 2;
		break;

	case 8:
		lColors = 256;
		BmpInfoHdr.biClrUsed = 256;
		break;

	case 24:
		lColors = 0;
		break;

	default:
		assert(FALSE);
	}
	
	lColorTableSize = lColors * sizeof(RGBQUAD);
	
	// Memory for 256 entries.

	lBISize = sizeof(BITMAPINFOHEADER) + lColorTableSize;
	lBitsSize = lWidthBytes * nZoomHeight;

	// Allocate memory.

	if ( ! ( hDib = GlobalAlloc(GMEM_MOVEABLE, lBISize + lBitsSize) ) )
		return FALSE;
	pBmpInfo = (LPBITMAPINFO) GlobalLock(hDib);
	memset(pBmpInfo, 0 , lBISize + lBitsSize);

	// Copy header to newly allocated memory.

	memcpy(pBmpInfo, &BmpInfoHdr, sizeof(BITMAPINFOHEADER));

	// Set the color table data for gray scale bitmaps.
    if ( lColors == 256 )
	{
		RGBQUAD* pRGB = (RGBQUAD*) (((LPBYTE)pBmpInfo) + sizeof(BITMAPINFOHEADER));
		for (LONG i=0; i<256; i++)
		{
			pRGB->rgbBlue = (BYTE)i;
    		pRGB->rgbGreen = (BYTE)i;
			pRGB->rgbRed = (BYTE)i;
			pRGB->rgbReserved = (BYTE)0;
			pRGB ++;
		}
	}
	
	// Set all members.
	// [(width*bitCount+31)/32]*4 对齐（不足32bit补齐），必须是4的倍数
	m_lWidthBytes = (m_pBMI->bmiHeader.biWidth * m_pBMI->bmiHeader.biBitCount + 31) / 32 * 4;

	BYTE* pZoomPixel = (BYTE*)(pBmpInfo) + lBISize;

	for(int y = 0; y < Height() - 1; y++)
		for(int x = 0; x < Width() - 1; x ++)
		{
			if( m_pBMI->bmiHeader.biBitCount == 8)
			{
				int nRowStart = y * ratio;
			    int nColStart = x * ratio;
			
				BYTE pixel1  = *((BYTE*)(m_pPixel + y * m_lWidthBytes + x));
				BYTE pixel2  = *((BYTE*)(m_pPixel + y * m_lWidthBytes + x + 1));
				BYTE pixel3  = *((BYTE*)(m_pPixel + (y + 1) * m_lWidthBytes + x));
				BYTE pixel4  = *((BYTE*)(m_pPixel + (y + 1) * m_lWidthBytes + x + 1));

				int nRatioSqr = ratio* ratio;
				for (int k = 0; k < ratio + 1; k++)
					for (int l = 0; l < ratio + 1; l++)
					{
						*(pZoomPixel + (nRowStart + k) * lWidthBytes + nColStart + l) =
							(BYTE)(( ((ratio - k) * (ratio - l) * pixel1) + 
									 ((ratio - k)* l * pixel2) +
									 ( k * (ratio - l) * pixel3) + 
								     ( k * l * pixel4)
								     ) / nRatioSqr);
					}
			}
			else if( m_pBMI->bmiHeader.biBitCount == 24)
			{
				int nRowStart = y * ratio;
			    int nColStart = x * ratio * 3;

				RGBTRIPLE pixel1  = *((RGBTRIPLE*)(m_pPixel + y * m_lWidthBytes + x *3));
				RGBTRIPLE pixel2  = *((RGBTRIPLE*)(m_pPixel + y * m_lWidthBytes + (x + 1) * 3));
				RGBTRIPLE pixel3  = *((RGBTRIPLE*)(m_pPixel + (y + 1) * m_lWidthBytes + x * 3));
				RGBTRIPLE pixel4  = *((RGBTRIPLE*)(m_pPixel + (y + 1) * m_lWidthBytes + (x + 1) *3));
				int nRatioSqr = ratio* ratio;
				
				for (int k = 0; k < ratio + 1; k++)
					for (int l = 0; l < ratio + 1; l++)
					{
						*(pZoomPixel + (nRowStart + k) * lWidthBytes + nColStart + l*3) =
							(BYTE)(( ((ratio - k) * (ratio - l) * pixel1.rgbtBlue) + 
									 ((ratio - k)* l * pixel2.rgbtBlue) +
									 ( k * (ratio - l) * pixel3.rgbtBlue) + 
								     ( k * l * pixel4.rgbtBlue)
								     ) / nRatioSqr);
					}
				for ( k = 0; k < ratio + 1; k++)
					for (int l = 0; l < ratio + 1; l++)
					{
						*(pZoomPixel + (nRowStart + k) * lWidthBytes + nColStart + l*3 + 1) =
							(BYTE)(( ((ratio - k) * (ratio - l) * pixel1.rgbtGreen) + 
									 ((ratio - k)* l * pixel2.rgbtGreen) +
									 ( k * (ratio - l) * pixel3.rgbtGreen) + 
								     ( k * l * pixel4.rgbtGreen)
								     ) / nRatioSqr);
					}
				for (k = 0; k < ratio + 1; k++)
					for (int l = 0; l < ratio + 1; l++)
					{
						*(pZoomPixel + (nRowStart + k) * lWidthBytes + nColStart + l*3 + 2) =
							(BYTE)(( ((ratio - k) * (ratio - l) * pixel1.rgbtRed) + 
									 ((ratio - k)* l * pixel2.rgbtRed) +
									 ( k * (ratio - l) * pixel3.rgbtRed) + 
								     ( k * l * pixel4.rgbtRed)
								     ) / nRatioSqr);
					}
			}
		}
	GlobalUnlock(hDib);
	return hDib;
}

HANDLE CImage::ConvertToGray()
{
	assert(m_pBMI->bmiHeader.biBitCount == 24);

	int width = Width();
	int height = Height();
	int bitCount = 8;
	
	HANDLE hDib;
	LPBITMAPINFO pBmpInfo;
	LONG lColors, lColorTableSize, lBitsSize, lBISize, lWidthBytes;

	BITMAPINFOHEADER BmpInfoHdr;

	BmpInfoHdr.biSize = sizeof(BITMAPINFOHEADER);
	BmpInfoHdr.biWidth = width;
	BmpInfoHdr.biHeight = height;
	BmpInfoHdr.biPlanes = 1;
	BmpInfoHdr.biBitCount = bitCount;
	BmpInfoHdr.biCompression = BI_RGB;
	BmpInfoHdr.biSizeImage = (width * bitCount + 31) / 32 * 4 * height;
	BmpInfoHdr.biXPelsPerMeter = 0;
	BmpInfoHdr.biYPelsPerMeter = 0;
	BmpInfoHdr.biClrUsed = 0;
	BmpInfoHdr.biClrImportant = 0;

	// The number of bytes for one scan line.

	lWidthBytes = (BmpInfoHdr.biWidth * BmpInfoHdr.biBitCount + 31) / 32 * 4;

	// Determine memory needed.

	switch (BmpInfoHdr.biBitCount)
	{
	case 1:
		lColors = 2;
		BmpInfoHdr.biClrUsed = 2;
		break;

	case 8:
		lColors = 256;
		BmpInfoHdr.biClrUsed = 256;
		break;

	case 24:
		lColors = 0;
		break;

	default:
		assert(FALSE);
	}
	
	lColorTableSize = lColors * sizeof(RGBQUAD);
	
	// Memory for 256 entries.

	lBISize = sizeof(BITMAPINFOHEADER) + lColorTableSize;
	lBitsSize = lWidthBytes * height;

	// Allocate memory.

	if ( ! ( hDib = GlobalAlloc(GMEM_MOVEABLE, lBISize + lBitsSize) ) )
		return FALSE;
	pBmpInfo = (LPBITMAPINFO) GlobalLock(hDib);

	// Copy header to newly allocated memory.

	memcpy(pBmpInfo, &BmpInfoHdr, sizeof(BITMAPINFOHEADER));

	// Set the color table data for gray scale bitmaps.

	if ( lColors == 256 )
	{
		RGBQUAD* pRGB = (RGBQUAD*) (((LPBYTE)pBmpInfo) + sizeof(BITMAPINFOHEADER));
		for (LONG i=0; i<256; i++)
		{
			pRGB->rgbBlue = (BYTE)i;
    		pRGB->rgbGreen = (BYTE)i;
			pRGB->rgbRed = (BYTE)i;
			pRGB->rgbReserved = (BYTE)0;
			pRGB ++;
		}
	}
	else if ( lColors == 2 )
	{
		RGBQUAD* pRGB = (RGBQUAD*) (((LPBYTE)pBmpInfo) + sizeof(BITMAPINFOHEADER));
		pRGB->rgbBlue = (BYTE) 0;
   		pRGB->rgbGreen = (BYTE) 0;
		pRGB->rgbRed = (BYTE) 0;
		pRGB->rgbReserved = (BYTE) 0;
		(pRGB + 1)->rgbBlue = (BYTE) 255;
   		(pRGB + 1)->rgbGreen = (BYTE) 255;
		(pRGB + 1)->rgbRed = (BYTE) 255;
		(pRGB + 1)->rgbReserved = (BYTE)0;
	}

	// Set all members.

	m_lWidthBytes = lWidthBytes * 3;
	BYTE* pGrayPixel = (BYTE*)(pBmpInfo) + lBISize;

	for(int y = 0; y < height; y++)
		for(int x = 0; x < width; x ++)
		{
			RGBTRIPLE rgbPixel = *((RGBTRIPLE*)(m_pPixel + y * m_lWidthBytes + 3 * x));
            // gray = r *0.114f + g*0.588f + b*0.298
			* (pGrayPixel + y * lWidthBytes + x) = 	(BYTE) 
				                (((int)rgbPixel.rgbtRed* 3735 + 
								 (int)rgbPixel.rgbtGreen * 19267 + 
								 (int)rgbPixel.rgbtBlue * 9765) 
								 /32767);
		}

	GlobalUnlock(hDib);
	return hDib;
}

HANDLE CImage::ExtractChannel(char cFlag)
{
	Lock();
	assert(m_pBMI->bmiHeader.biBitCount == 24);

	int width = Width();
	int height = Height();
	int bitCount = 8;
	
	HANDLE hDib;
	LPBITMAPINFO pBmpInfo;
	LONG lColors, lColorTableSize, lBitsSize, lBISize, lWidthBytes;

	BITMAPINFOHEADER BmpInfoHdr;

	BmpInfoHdr.biSize = sizeof(BITMAPINFOHEADER);
	BmpInfoHdr.biWidth = width;
	BmpInfoHdr.biHeight = height;
	BmpInfoHdr.biPlanes = 1;
	BmpInfoHdr.biBitCount = bitCount;
	BmpInfoHdr.biCompression = BI_RGB;
	BmpInfoHdr.biSizeImage = (width * bitCount + 31) / 32 * 4 * height;
	BmpInfoHdr.biXPelsPerMeter = 0;
	BmpInfoHdr.biYPelsPerMeter = 0;
	BmpInfoHdr.biClrUsed = 0;
	BmpInfoHdr.biClrImportant = 0;

	// The number of bytes for one scan line.

	lWidthBytes = (BmpInfoHdr.biWidth * BmpInfoHdr.biBitCount + 31) / 32 * 4;

	// Determine memory needed.

	switch (BmpInfoHdr.biBitCount)
	{
	case 1:
		lColors = 2;
		BmpInfoHdr.biClrUsed = 2;
		break;

	case 8:
		lColors = 256;
		BmpInfoHdr.biClrUsed = 256;
		break;

	case 24:
		lColors = 0;
		break;

	default:
		assert(FALSE);
	}
	
	lColorTableSize = lColors * sizeof(RGBQUAD);
	
	// Memory for 256 entries.

	lBISize = sizeof(BITMAPINFOHEADER) + lColorTableSize;
	lBitsSize = lWidthBytes * height;

	// Allocate memory.

	if ( ! ( hDib = GlobalAlloc(GMEM_MOVEABLE, lBISize + lBitsSize) ) )
		return FALSE;
	pBmpInfo = (LPBITMAPINFO) GlobalLock(hDib);

	// Copy header to newly allocated memory.

	memcpy(pBmpInfo, &BmpInfoHdr, sizeof(BITMAPINFOHEADER));

	// Set the color table data for gray scale bitmaps.

	if ( lColors == 256 )
	{
		RGBQUAD* pRGB = (RGBQUAD*) (((LPBYTE)pBmpInfo) + sizeof(BITMAPINFOHEADER));
		for (LONG i=0; i<256; i++)
		{
			pRGB->rgbBlue = (BYTE)i;
    		pRGB->rgbGreen = (BYTE)i;
			pRGB->rgbRed = (BYTE)i;
			pRGB->rgbReserved = (BYTE)0;
			pRGB ++;
		}
	}
	else if ( lColors == 2 )
	{
		RGBQUAD* pRGB = (RGBQUAD*) (((LPBYTE)pBmpInfo) + sizeof(BITMAPINFOHEADER));
		pRGB->rgbBlue = (BYTE) 0;
   		pRGB->rgbGreen = (BYTE) 0;
		pRGB->rgbRed = (BYTE) 0;
		pRGB->rgbReserved = (BYTE) 0;
		(pRGB + 1)->rgbBlue = (BYTE) 255;
   		(pRGB + 1)->rgbGreen = (BYTE) 255;
		(pRGB + 1)->rgbRed = (BYTE) 255;
		(pRGB + 1)->rgbReserved = (BYTE)0;
	}

	

	// Set all members.

	m_lWidthBytes = lWidthBytes * 3;

	BYTE* pGrayPixel = (BYTE*)(pBmpInfo) + lBISize;

	for(int y = 0; y < height; y++)
		for(int x = 0; x < width; x ++)
		{
			RGBTRIPLE rgbPixel = *((RGBTRIPLE*)(m_pPixel + y * m_lWidthBytes + 3 * x));
            switch (cFlag)
			{
			case 'R':
				* (pGrayPixel + y * lWidthBytes + x) = (BYTE)(rgbPixel.rgbtRed);
				break;
            case 'G':
                * (pGrayPixel + y * lWidthBytes + x) = (BYTE)(rgbPixel.rgbtGreen);
				break;
			case 'B':
				* (pGrayPixel + y * lWidthBytes + x) = (BYTE)(rgbPixel.rgbtBlue);
				break;
			default:
                MessageBox(NULL, "Invalidate Channel.", "Channel specification is invalidate", MB_OK);
         
			}
		}
	GlobalUnlock(hDib);
	Unlock();
	return hDib;

}
/////////////////////////////////////////////////////////////////////////////

BOOL CImage::Create(LONG width, LONG height, WORD bitCount)
{
	HANDLE hDib;
	LPBITMAPINFO pBmpInfo;
	LONG lColors, lColorTableSize, lBitsSize, lBISize, lWidthBytes;

	BITMAPINFOHEADER BmpInfoHdr;

	BmpInfoHdr.biSize = sizeof(BITMAPINFOHEADER);
	BmpInfoHdr.biWidth = width;
	BmpInfoHdr.biHeight = height;
	BmpInfoHdr.biPlanes = 1;
	BmpInfoHdr.biBitCount = bitCount;
	BmpInfoHdr.biCompression = BI_RGB;
	BmpInfoHdr.biSizeImage = (width * bitCount + 31) / 32 * 4 * height;
	BmpInfoHdr.biXPelsPerMeter = 0;
	BmpInfoHdr.biYPelsPerMeter = 0;
	BmpInfoHdr.biClrUsed = 0;
	BmpInfoHdr.biClrImportant = 0;

	// The number of bytes for one scan line.

	lWidthBytes = (BmpInfoHdr.biWidth * BmpInfoHdr.biBitCount + 31) / 32 * 4;

	// Determine memory needed.

	switch (BmpInfoHdr.biBitCount)
	{
	case 1:
		lColors = 2;
		BmpInfoHdr.biClrUsed = 2;
		break;

	case 8:
		lColors = 256;
		BmpInfoHdr.biClrUsed = 256;
		break;
	case 16:
	case 24:
		lColors = 0;
		break;

	default:
		assert(FALSE);
	}
	
	lColorTableSize = lColors * sizeof(RGBQUAD);
	
	// Memory for 256 entries.

	lBISize = sizeof(BITMAPINFOHEADER) + lColorTableSize;
	lBitsSize = lWidthBytes * height;

	// Allocate memory.

	if ( ! ( hDib = GlobalAlloc(GMEM_MOVEABLE, lBISize + lBitsSize) ) )
//	if ( ! ( hDib = GlobalAlloc(GHND, lBISize + lBitsSize) ) )
		return FALSE;
	pBmpInfo = (LPBITMAPINFO) GlobalLock(hDib);

	// Copy header to newly allocated memory.

	memcpy(pBmpInfo, &BmpInfoHdr, sizeof(BITMAPINFOHEADER));

	// Set the color table data for gray scale bitmaps.

	if ( lColors == 256 )
	{
		RGBQUAD* pRGB = (RGBQUAD*) (((LPBYTE)pBmpInfo) + sizeof(BITMAPINFOHEADER));
		for (LONG i=0; i<256; i++)
		{
			pRGB->rgbBlue = (BYTE)i;
    		pRGB->rgbGreen = (BYTE)i;
			pRGB->rgbRed = (BYTE)i;
			pRGB->rgbReserved = (BYTE)0;
			pRGB ++;
		}
	}
	else if ( lColors == 2 )
	{
		RGBQUAD* pRGB = (RGBQUAD*) (((LPBYTE)pBmpInfo) + sizeof(BITMAPINFOHEADER));
		pRGB->rgbBlue = (BYTE) 0;
   		pRGB->rgbGreen = (BYTE) 0;
		pRGB->rgbRed = (BYTE) 0;
		pRGB->rgbReserved = (BYTE) 0;
		(pRGB + 1)->rgbBlue = (BYTE) 255;
   		(pRGB + 1)->rgbGreen = (BYTE) 255;
		(pRGB + 1)->rgbRed = (BYTE) 255;
		(pRGB + 1)->rgbReserved = (BYTE)0;
	}

	GlobalUnlock(hDib);

	// Set all members.

	Release();

	m_hDib = hDib;
	m_lWidthBytes = lWidthBytes;

	m_pBMI = NULL;
	m_pPixel = NULL;

	m_bIsAttached = FALSE;

	return TRUE;
}


VOID CImage::Release()
{
	if ( m_bIsAttached )
		Detach();
	else
		if ( m_hDib )
			GlobalFree(m_hDib);	// Warning: Return value ignored.

	m_hDib = NULL;
}


/////////////////////////////////////////////////////////////////////////////

BOOL CImage::Attach(HANDLE hDib)
{
	LPBITMAPINFO pBMI;

	if ( ! ( pBMI = (LPBITMAPINFO) GlobalLock(hDib) ) )
		return FALSE;

	// 1/8/24 bits images only

//	assert(pBMI->bmiHeader.biBitCount == 1 || pBMI->bmiHeader.biBitCount == 8 || pBMI->bmiHeader.biBitCount == 24);

	pBMI->bmiHeader.biSizeImage = 
		(pBMI->bmiHeader.biWidth * pBMI->bmiHeader.biBitCount + 31) / 32 * 4 * pBMI->bmiHeader.biHeight;

	// The number of bytes for one scan line.
	m_lWidthBytes = (pBMI->bmiHeader.biWidth * pBMI->bmiHeader.biBitCount + 31) / 32 * 4;

	GlobalUnlock(hDib);

	Release();

	m_hDib = hDib;
	m_pBMI = NULL;
	m_pPixel = NULL;
	m_bIsAttached = TRUE;

	return TRUE;
}

HANDLE CImage::Detach()
{
	HANDLE hDib = m_hDib;
	m_hDib = NULL;
	m_pBMI = NULL;
	m_pPixel = NULL;
	m_bIsAttached = FALSE;

	return hDib;
}


/////////////////////////////////////////////////////////////////////////////

// Draw the DIB to the specified device context at position (x,y).

VOID CImage::Draw(HDC hDC, LONG x, LONG y)
{
	assert(m_pBMI);

	::SetDIBitsToDevice(
		hDC,						// DC handle
		x,							// Destination x
		y,							// Destination y
		m_pBMI->bmiHeader.biWidth,	// Source width
		m_pBMI->bmiHeader.biHeight,	// Source height
		0,							// Source x
		0,							// Source y
		0,							// first scan line in array 
		m_pBMI->bmiHeader.biHeight,	// number of scan lines 
		m_pPixel,					// Pointer to bits
		m_pBMI,						// Pointer to BITMAPINFO
		DIB_RGB_COLORS				// Options
	);	
}


// Stretch the DIB to the specified device context.

VOID CImage::Stretch(HDC hDC, LONG x, LONG y, LONG width, LONG height)
{
	assert(m_pBMI);

	::StretchDIBits(
		hDC,						// DC handle
		x,							// Destination x
		y,							// Destination y
		width,						// Destination width
		height,						// Destination height
		0,							// Source x
		0,							// Source y
		m_pBMI->bmiHeader.biWidth,	// Source width
		m_pBMI->bmiHeader.biHeight,	// Source height
		m_pPixel,					// Pointer to bits
		m_pBMI,						// Pointerer to BITMAPINFO
		DIB_RGB_COLORS,				// Options,
		SRCCOPY);					// Raster operation code (ROP)
}


/////////////////////////////////////////////////////////////////////////////

VOID CImage::Copy(CImage* pDestDIB)
{
	assert(pDestDIB);
	assert(pDestDIB->Width() == Width());
	assert(pDestDIB->Height() == Height());

	memcpy(pDestDIB->Pixel(), Pixel(), BMI()->bmiHeader.biSizeImage);
}


/////////////////////////////////////////////////////////////////////////////

VOID CImage::Clear(BYTE byInit)
{
	memset( m_pPixel, byInit, Size() );
}


/////////////////////////////////////////////////////////////////////////////

BOOL CImage::Gaussian(CImage* pDest)
{
	if ( m_pBMI->bmiHeader.biBitCount == 24 )
		return GaussianRGB(pDest);

	LONG lRow = Height();
	LONG lCol = WidthBytes();

	HANDLE hBuf = GlobalAlloc(GMEM_MOVEABLE, lRow * lCol);
	if ( ! hBuf )
		return FALSE;
	LPBYTE pBuf = (LPBYTE) GlobalLock(hBuf);

	BYTE* pSrc = Pixel();
	BYTE* pDst = pDest->Pixel();
	BYTE *pB, *pS, *pD;
	LONG x, y;

	for (y=lRow - 1; y>=0; y--)
	{
		pB = pBuf + y * lCol;
		pS = pSrc + y * lCol;
		*pB = ( *pS + *pS + *(pS + 1) ) / 3;
		pB ++;  pS ++;
		for (x = lCol - 2; x>=1; x--)
		{
			*pB = ( *(pS-1) + *pS + *pS + *(pS+1) ) >> 2;
			pB ++;  pS ++;
		}
		*pB = ( *(pS-1) + *pS + *pS ) / 3;
	}

	pB = pBuf;
	pD = pDst;
	for (x=lCol-1; x>=0; x--)
	{
		*pD = ( *pB + *pB + *(pB + lCol) ) / 3;
		pB ++;  pD ++;
	}

	pB = pBuf + (lRow - 1) * lCol;
	pD = pDst + (lRow - 1) * lCol;
	for (x=lCol - 1; x>=0; x--)
	{
		*pD = ( *(pB - lCol) + *pB + *pB ) / 3;
		pB ++;  pD ++;
	}

	for (y=lRow - 2; y>=1; y--)
	{
		pB = pBuf + y*lCol;
		pD = pDst + y*lCol;
		for (x=lCol - 1; x>=0; x--)
		{
			*pD = ( *(pB - lCol) + *pB + *pB + *(pB + lCol) ) >> 2;
			pB ++;  pD ++;
		}
	}

	GlobalUnlock(hBuf);
	GlobalFree(hBuf);

	return TRUE;
}

BOOL CImage::GaussianRGB(CImage* pDest)
{
	assert(m_pBMI->bmiHeader.biBitCount == 24);

	LONG cx = Width();
	LONG cy = Height();

	LONG x, y;

	HANDLE hBuf = GlobalAlloc(GMEM_MOVEABLE, cx * cy * sizeof(RGBTRIPLE));
	if ( ! hBuf )
		return FALSE;
	RGBTRIPLE* pBuf = (RGBTRIPLE*) GlobalLock(hBuf);

	RGBTRIPLE rgbt0, rgbt1, rgbt2, rgbt;

	for (y=cy - 1; y>=0; y--)
	{
		rgbt0 = GetRGB(0, y);
		rgbt1 = GetRGB(1, y);
		pBuf[y * cx + 0].rgbtRed = (2 * rgbt0.rgbtRed + rgbt1.rgbtRed) / 3;
		pBuf[y * cx + 0].rgbtGreen = (2 * rgbt0.rgbtGreen + rgbt1.rgbtGreen) / 3;
		pBuf[y * cx + 0].rgbtBlue = (2 * rgbt0.rgbtBlue + rgbt1.rgbtBlue) / 3;

		for (x=cx - 2; x>=1; x--)
		{
			rgbt0 = GetRGB(x, y);
			rgbt1 = GetRGB(x - 1, y);
			rgbt2 = GetRGB(x + 1, y);
			pBuf[y * cx + x].rgbtBlue = (((rgbt0.rgbtBlue << 1) + rgbt1.rgbtBlue + rgbt2.rgbtBlue) >> 2);
			pBuf[y * cx + x].rgbtGreen = (((rgbt0.rgbtGreen << 1) + rgbt1.rgbtGreen + rgbt2.rgbtGreen) >> 2);
			pBuf[y * cx + x].rgbtRed = (((rgbt0.rgbtRed << 1) + rgbt1.rgbtRed + rgbt2.rgbtRed) >> 2);
		}

		rgbt0 = GetRGB(cx - 1, y);
		rgbt1 = GetRGB(cx - 2, y);
		pBuf[y * cx + cx - 1].rgbtRed = (2 * rgbt0.rgbtRed + rgbt1.rgbtRed) / 3;
		pBuf[y * cx + cx - 1].rgbtGreen = (2 * rgbt0.rgbtGreen + rgbt1.rgbtGreen) / 3;
		pBuf[y * cx + cx - 1].rgbtBlue = (2 * rgbt0.rgbtBlue + rgbt1.rgbtBlue) / 3;
	}
	
	for (x=cx - 1; x>=0; x--)
	{
		rgbt0 = pBuf[0 * cy + x];
		rgbt1 = pBuf[1 * cy + x];
		rgbt.rgbtRed = (2 * rgbt0.rgbtRed + rgbt1.rgbtRed) / 3;
		rgbt.rgbtGreen = (2 * rgbt0.rgbtGreen + rgbt1.rgbtGreen) / 3;
		rgbt.rgbtBlue = (2 * rgbt0.rgbtBlue + rgbt1.rgbtBlue) / 3;
		pDest->PutRGB(x, 0, rgbt);

		for (y=cy - 2; y>=1; y--)
		{
			rgbt0 = pBuf[y * cx + x];
			rgbt1 = pBuf[(y - 1) * cx + x];
			rgbt2 = pBuf[(y + 1) * cx + x];
			rgbt.rgbtRed = (((rgbt0.rgbtRed << 1) + rgbt1.rgbtRed + rgbt2.rgbtRed) >> 2);
			rgbt.rgbtGreen = (((rgbt0.rgbtGreen << 1) + rgbt1.rgbtGreen + rgbt2.rgbtGreen) >> 2);
			rgbt.rgbtBlue = (((rgbt0.rgbtBlue << 1) + rgbt1.rgbtBlue + rgbt2.rgbtBlue) >> 2);
			pDest->PutRGB(x, y, rgbt);
		}

		rgbt0 = pBuf[(cy - 1) * cx + x];
		rgbt1 = pBuf[(cy - 2) * cx + x];
		rgbt.rgbtRed = (2 * rgbt0.rgbtRed + rgbt1.rgbtRed) / 3;
		rgbt.rgbtGreen = (2 * rgbt0.rgbtGreen + rgbt1.rgbtGreen) / 3;
		rgbt.rgbtBlue = (2 * rgbt0.rgbtBlue + rgbt1.rgbtBlue) / 3;
		pDest->PutRGB(x, cy - 1, rgbt);
	}

	GlobalUnlock(hBuf);
	GlobalFree(hBuf);

	return TRUE;
}

BOOL CImage::LoadFromCurPos(FILE* pfile)
{
	HANDLE hBmpInfo = NULL;
	BITMAPINFO* pBmpInfo;
	BYTE* pBits = NULL;
	LONG lColors, lColorTableSize, lBitsSize, lBISize;

	// Read the file header.  This will tell us the file size and
	// where the data bits start in the file.
	
	BITMAPFILEHEADER BmpFileHdr;
	LONG lBytes;
	lBytes = fread(&BmpFileHdr, 1, sizeof(BITMAPFILEHEADER), pfile);
	if (lBytes != sizeof(BITMAPFILEHEADER))
		goto error;

	// Do we have "BM" at the start indicating a bitmap file?

	if (BmpFileHdr.bfType != 0x4D42)
		goto error;

	// Assume for now that the file is a Wondows DIB.
	// Read the BITMAPINFOHEADER.

	BITMAPINFOHEADER BmpInfoHdr;
	lBytes = fread(&BmpInfoHdr, 1, sizeof(BITMAPINFOHEADER), pfile);
	if (lBytes != sizeof(BITMAPINFOHEADER))
		goto error;

	// Fill the biSizeImage field.

	if ( BmpInfoHdr.biSizeImage == 0 )
		BmpInfoHdr.biSizeImage = ( BmpInfoHdr.biWidth * BmpInfoHdr.biBitCount + 31 ) / 32 * 4 * BmpInfoHdr.biHeight;

	// The number of bytes for one scan line.

	m_lWidthBytes = (BmpInfoHdr.biWidth * BmpInfoHdr.biBitCount + 31) / 32 * 4;

	// Was it a real Windows DIB file?

	if (BmpInfoHdr.biSize != sizeof(BITMAPINFOHEADER))
		goto error;

	// Determine memory needed.

	lColors = GetColorEntries((LPBITMAPINFO)&BmpInfoHdr);
	lColorTableSize = lColors * sizeof(RGBQUAD);
	
	// Allocate memory for color table entries.
	
	lBISize = sizeof(BITMAPINFOHEADER) + lColorTableSize;
	lBitsSize = BmpFileHdr.bfSize - BmpFileHdr.bfOffBits;

	// Allocate memory.

	if ( ! ( hBmpInfo = GlobalAlloc(GMEM_MOVEABLE, lBISize + lBitsSize) ) )
		goto error;
	pBmpInfo = (LPBITMAPINFO) GlobalLock(hBmpInfo);

	// Copy header to newly allocated memory.
	
	memcpy(pBmpInfo, &BmpInfoHdr, sizeof(BITMAPINFOHEADER));

	// Read the color table data from the file.
	
	if ( lColorTableSize > 0 )
	{
		lBytes = fread(((LPBYTE)pBmpInfo) + sizeof(BITMAPINFOHEADER), 1, lColorTableSize, pfile);
		if (lBytes != lColorTableSize)
			goto error;
	}

	// Set bits

	pBits = (LPBYTE) pBmpInfo + lBISize;

	// Read data.
	
	lBytes = fread(pBits, 1, lBitsSize, pfile);
	if (lBytes != lBitsSize)
		goto error;

	// Fill the biSizeImage field.
	
	pBmpInfo->bmiHeader.biSizeImage = lBitsSize;

	// We reach here only if if no errors occured.

	GlobalUnlock(hBmpInfo);

	Release();

	m_hDib = hBmpInfo;
	m_bIsAttached = FALSE;

	fclose(pfile);

	return TRUE;
	
	// Error manipulations.

error:

	if ( hBmpInfo)
	{
		GlobalUnlock(hBmpInfo);
		GlobalFree(hBmpInfo);
	}

	fclose(pfile);
}
/////////////////////////////////////////////////////////////////////////////

BOOL CImage::Load(LPCTSTR lpszFilePath)
{
	// Open as a bitmap 

	FILE* f = fopen(lpszFilePath, "rb");
	
	if ( ! f )
	{
		printf("Bitmap File Load Failed!");
		return FALSE;
	}

	HANDLE hBmpInfo = NULL;
	BITMAPINFO* pBmpInfo;
	BYTE* pBits = NULL;
	LONG lColors, lColorTableSize, lBitsSize, lBISize;

	// Get the current file position.

	LONG lFileStart = ftell(f);

	// Read the file header.  This will tell us the file size and
	// where the data bits start in the file.
	
	BITMAPFILEHEADER BmpFileHdr;
	LONG lBytes;
	lBytes = fread(&BmpFileHdr, 1, sizeof(BITMAPFILEHEADER), f);
	if (lBytes != sizeof(BITMAPFILEHEADER))
		goto error;

	// Do we have "BM" at the start indicating a bitmap file?

	if (BmpFileHdr.bfType != 0x4D42)
		goto error;

	// Assume for now that the file is a Wondows DIB.
	// Read the BITMAPINFOHEADER.

	BITMAPINFOHEADER BmpInfoHdr;
	lBytes = fread(&BmpInfoHdr, 1, sizeof(BITMAPINFOHEADER), f);
	if (lBytes != sizeof(BITMAPINFOHEADER))
		goto error;

	// Fill the biSizeImage field.
	if(BmpInfoHdr.biBitCount != 24)
		return FALSE;

	m_BMPwidth = BmpInfoHdr.biWidth;
    m_BMPheigth = BmpInfoHdr.biHeight;

	if ( BmpInfoHdr.biSizeImage == 0 )//////////////////////////////////////////////////////
		BmpInfoHdr.biSizeImage = ( BmpInfoHdr.biWidth * BmpInfoHdr.biBitCount + 31 ) / 32 * 4 * BmpInfoHdr.biHeight;

	// The number of bytes for one scan line.

	m_lWidthBytes = (BmpInfoHdr.biWidth * BmpInfoHdr.biBitCount + 31) / 32 * 4;

	// Was it a real Windows DIB file?

	if (BmpInfoHdr.biSize != sizeof(BITMAPINFOHEADER))
		goto error;

	// Determine memory needed.

	lColors = GetColorEntries((LPBITMAPINFO)&BmpInfoHdr);
	lColorTableSize = lColors * sizeof(RGBQUAD);
	
	// Allocate memory for color table entries.
	
	lBISize = sizeof(BITMAPINFOHEADER) + lColorTableSize;
	lBitsSize = BmpFileHdr.bfSize - BmpFileHdr.bfOffBits;

	// Allocate memory.

	if ( ! ( hBmpInfo = GlobalAlloc(GMEM_MOVEABLE, lBISize + lBitsSize) ) )
		goto error;
	pBmpInfo = (LPBITMAPINFO) GlobalLock(hBmpInfo);

	// Copy header to newly allocated memory.
	
	memcpy(pBmpInfo, &BmpInfoHdr, sizeof(BITMAPINFOHEADER));

	// Read the color table data from the file.
	
	if ( lColorTableSize > 0 )
	{
		lBytes = fread(((LPBYTE)pBmpInfo) + sizeof(BITMAPINFOHEADER), 1, lColorTableSize, f);
		if (lBytes != lColorTableSize)
			goto error;
	}

	// Set bits

	pBits = (LPBYTE) pBmpInfo + lBISize;

	// Move file pointer to file location of bits.
	
	fseek(f, lFileStart + BmpFileHdr.bfOffBits, SEEK_SET);

	// Read data.
	
	lBytes = fread(pBits, 1, lBitsSize, f);
	if (lBytes != lBitsSize)
		goto error;

	// Fill the biSizeImage field.
	
	pBmpInfo->bmiHeader.biSizeImage = lBitsSize;

	// We reach here only if if no errors occured.

	GlobalUnlock(hBmpInfo);

	Release();

	m_hDib = hBmpInfo;
	m_bIsAttached = FALSE;

	fclose(f);

	return TRUE;
	
	// Error manipulations.

error:

	if ( hBmpInfo)
	{
		GlobalUnlock(hBmpInfo);
		GlobalFree(hBmpInfo);
	}

	fclose(f);

	return FALSE;
}


// Save DIB to a file.
BOOL CImage::AddtoFileEnd(FILE *pfile)
{
	BITMAPFILEHEADER BmpFileHdr;
	BITMAPINFOHEADER* pBmpInfoHdr = (LPBITMAPINFOHEADER) m_pBMI;
	LONG lColors, lColorTableSize, lBitsSize, lBISize;

	// Determine size.

	lColors = GetColorEntries(m_pBMI);
	lColorTableSize = lColors * sizeof(RGBQUAD);
	lBISize = sizeof(BITMAPINFOHEADER) + lColorTableSize;
	DWORD dwWidthBytes = (pBmpInfoHdr->biWidth*pBmpInfoHdr->biBitCount + 31) / 32 * 4;
	lBitsSize = dwWidthBytes * pBmpInfoHdr->biHeight;

	// Write the file header.

	BmpFileHdr.bfType = 0x4D42;		// 'BM'
	BmpFileHdr.bfSize = sizeof(BITMAPFILEHEADER) + lBISize + lBitsSize;
	BmpFileHdr.bfReserved1 = 0;
	BmpFileHdr.bfReserved2 = 0;
	BmpFileHdr.bfOffBits = sizeof(BITMAPFILEHEADER) + lBISize;

	fwrite(&BmpFileHdr, sizeof(BITMAPFILEHEADER), 1, pfile);

	// Write the BITMAPINFOHEADER and color table.

	fwrite(m_pBMI, lBISize, 1, pfile);

	// Write bits data.

	fwrite(m_pPixel,lBitsSize, 1, pfile);

	return TRUE;
}

BOOL CImage::Save(LPCTSTR lpszFilePath)
{
	FILE* f = fopen(lpszFilePath, "wb+");

	Lock();

	BITMAPFILEHEADER BmpFileHdr;
	BITMAPINFOHEADER* pBmpInfoHdr = (LPBITMAPINFOHEADER) m_pBMI;
	LONG lColors, lColorTableSize, lBitsSize, lBISize;

	// Determine size.

	lColors = GetColorEntries(m_pBMI);
	lColorTableSize = lColors * sizeof(RGBQUAD);
	lBISize = sizeof(BITMAPINFOHEADER) + lColorTableSize;
	DWORD dwWidthBytes = (pBmpInfoHdr->biWidth*pBmpInfoHdr->biBitCount + 31) / 32 * 4;
	lBitsSize = dwWidthBytes * pBmpInfoHdr->biHeight;

	// Write the file header.

	BmpFileHdr.bfType = 0x4D42;		// 'BM'
	BmpFileHdr.bfSize = sizeof(BITMAPFILEHEADER) + lBISize + lBitsSize;
	BmpFileHdr.bfReserved1 = 0;
	BmpFileHdr.bfReserved2 = 0;
	BmpFileHdr.bfOffBits = sizeof(BITMAPFILEHEADER) + lBISize;

	fwrite(&BmpFileHdr, 1, sizeof(BITMAPFILEHEADER), f);

	// Write the BITMAPINFOHEADER and color table.

	fwrite(m_pBMI, 1, lBISize, f);

	// Write bits data.

	fwrite(m_pPixel, 1, lBitsSize, f);

	Unlock();

	fclose(f);

	return TRUE;
}

double CImage::GetPSNRColor(CImage* pRecon,int iBW,int iBH) // 24bits bmp
{
	CImage* pSrc = this;

	pSrc->Lock();
	pRecon->Lock();

	int iWidth  = pSrc->Width();
	int iHeight = pSrc->Height();
	
	HANDLE hSrc = pSrc->ConvertToGray();
	HANDLE hRecon = pRecon->ConvertToGray();

	pSrc->Unlock();
	pRecon->Unlock();


	CImage imgSrc,imgRecon;
	imgSrc.Attach(hSrc);
	imgRecon.Attach(hRecon);

	double psnr = imgSrc.GetPSNRGray(&imgRecon,iBW,iBH);

	imgSrc.Detach();
	imgRecon.Detach();
	GlobalFree(hSrc);
	GlobalFree(hRecon);

	return psnr;
}

double CImage::GetPSNRGray(CImage* pRecon,int iBW,int iBH) 
{
	CImage* pSrc = this;
	
	pSrc->Lock();
	pRecon->Lock();

	int iWidth = pSrc->Width();
	int iHeight = pSrc->Height();

	assert(pSrc->BMI()->bmiHeader.biBitCount == 8);

	int iLeft = iBW;
    int iRight = iWidth - 1 - iBW;
	int iBottom = iBH;
	int iTop = iHeight - 1 - iBH;

	//iLeft = (iLeft < iWidth):iLeft ? 0;
	//iRight = (iRight >=0) : iRight ? 0;
	//iBottom = (iBottom < iHeight) : iBottom ? 0;
	//iTop = (iTop >=0) : iTop ? 0;

	int w = iRight - iLeft +1;
	int h = iTop - iBottom +1;
	int area = w*h;
	
	int iSum  = 0;
	for(int cy = iBottom; cy <= iTop; cy++)
	{
		for(int cx  = iLeft; cx <= iRight; cx++)
		{
			int iErr =  pSrc->Get(cx,cy) - pRecon->Get(cx,cy);
			iSum  = iSum + iErr * iErr;
		}
	}
	double mse = (double)iSum/area;
	double psnr = 10* log10((255*255)/mse);

	pSrc->Unlock();
	pRecon->Unlock();
	
	return psnr;
}

BOOL CImage::GetYIQ(float* pYData,float*pIData,float* pQData)
{
	HANDLE hDib = m_hDib;
	CImage img;
	
	img.Attach(hDib);
	img.Lock();
	
	int iWidth = img.Width();
	int iHeight = img.Height();
	
	assert(img.BMI()->bmiHeader.biBitCount == 24);
	
	for(int cy = 0; cy < iHeight; cy++)
	{
		for(int cx = 0; cx < iWidth; cx++)
		{
			RGBTRIPLE rgb = img.GetRGB(cx,cy);
			*pYData = ( 299 * rgb.rgbtRed + 587 * rgb.rgbtGreen + 114 * rgb.rgbtBlue)/1000;
			*pIData = ( 596 * rgb.rgbtRed - 275 * rgb.rgbtGreen - 321 * rgb.rgbtBlue)/1000;
			*pQData = ( 212 * rgb.rgbtRed - 523 * rgb.rgbtGreen + 311 * rgb.rgbtBlue)/1000;
		
			pYData++;
			pIData++;
			pQData++;
		}
	}
	
	img.Unlock();
	img.Detach();
	return TRUE;
}

CImage* CImage::CombineYIQ(float* pYData,float* pIData, float* pQData)
{
	int iWidth = this->Width();
	int iHeight = this->Height();

	CImage* pImg = new CImage();
	pImg->Create(iWidth,iHeight,24);
	
	pImg->Lock();
	for(int cy = 0; cy < iHeight; cy++)
	for(int cx = 0; cx < iWidth; cx++)
	{
        int Y = (int)(*pYData);
		int I = (int)(*pIData);
		int Q = (int)(*pQData);
	    		
		int r = ( 1000 * Y + 956 * I + 621 * Q)/1000;
	    int g = ( 1000 * Y - 272 * I - 647 * Q)/1000;
        int b = ( 1000 * Y - 1107 * I + 1704 * Q)/1000;
		
		r = (r < 0) ? 0 : r;
        r = (r > 255) ? 255 : r;
        
		g = (g < 0) ? 0 : g;
        g = (g > 255) ? 255 : g;
		b = (b < 0) ? 0 : b;
        b = (b > 255) ? 255 : b;

		RGBTRIPLE rgb;
		rgb.rgbtRed = r;
		rgb.rgbtGreen = g;
		rgb.rgbtBlue = b;
        
		pImg->PutRGB(cx,cy,rgb);
   
		pYData++;
		pIData++;
		pQData++;
	}
	pImg->Unlock();
	return pImg;
}

CImage* CImage::TransToYIQ()
{
	Lock();
	int iWidth = this->Width();
	int iHeight = this->Height();
    
	float* pYData = new float[iWidth*iHeight];
	float* pIData = new float[iWidth*iHeight];
	float* pQData = new float[iWidth*iHeight];
    
	GetYIQ(pYData,pIData,pQData);
	CImage* pResult = CombineYIQ(pYData,pIData,pQData);

	delete[] pYData;
	delete[] pIData;
	delete[] pQData;
    Unlock();

	return pResult;
}

BOOL CImage::GetPixel(int x, int y, COLORREF &value)
{
	if(m_hDib == NULL)
		return FALSE;

	if(GetBitsPerPixel() != 24)
		return FALSE;

	if (x >= Width() || x < 0 || y >= Height() || y < 0)
		return FALSE;

	BYTE* pData = m_pPixel + (Height() - 1 - y)*m_lWidthBytes + x*3;
	value = RGB(pData[2], pData[1], pData[0]);

	return TRUE;
}

int CImage::GetBitsPerPixel() const
{
	int nBitCount;

	if (m_hDib == NULL)
		nBitCount = 0;
	else
	{
  		LPSTR lpDib = (LPSTR) ::GlobalLock(m_hDib);

		//  Calculate the number of bits per pixel for the DIB.
		if (IS_WIN30_DIB(lpDib))
			nBitCount = ((LPBITMAPINFOHEADER)lpDib)->biBitCount;
		else
			nBitCount = ((LPBITMAPCOREHEADER)lpDib)->bcBitCount;

	 	::GlobalUnlock((HGLOBAL) m_hDib);
	}

	return nBitCount;
}

BOOL CImage::SetPixel(int x, int y, const COLORREF &value)
{
	int bitCount = GetBitsPerPixel();
	BYTE* pData;

	if (m_hDib == NULL)
		return FALSE;

	if(bitCount != 24)
		return FALSE;
	//Currently only support pixel access for 24 bit images
	//May support other formats in the future
	//check the position
	if(!(x < Width() && x >= 0 && y < Height() && y >= 0))printf(" x , y has gone beyond width/height./n");
	
	pData = m_pPixel + (Height()-1-y)*m_lWidthBytes + x*3;  //THIS NEEDS TO BE TESTED
		
	pData[0] = GetBValue(value);
	pData[1] = GetGValue(value);
	pData[2] = GetRValue(value);

	return TRUE;
}

void CImage::RGBtoHSL(COLORREF rgb, double *H, double *S, double *L)
{
	double delta;
	double r = (double)GetRValue(rgb)/255;
	double g = (double)GetGValue(rgb)/255;
	double b = (double)GetBValue(rgb)/255;
	double cmax = max(r,max(g,b));
	double cmin = min(r,min(g,b));

	*L = (cmax+cmin)/2.0;
	if (cmax == cmin) 
	{
		*S = 0;
		*H = 0; // it's really undefined
	} 
	else 
	{
		if (*L < 0.5) 
			*S = (cmax-cmin)/(cmax+cmin);
		else
			*S = (cmax-cmin)/(2.0-cmax-cmin);
		delta = cmax - cmin;
		if (r==cmax)
			*H = (g-b)/delta;
		else if (g==cmax)
			*H = 2.0 +(b-r)/delta;
		else
			*H = 4.0+(r-g)/delta;
		*H /= 6.0;
		if (*H < 0.0)
			*H += 1;
	}
}

double CImage::HuetoRGB(double m1, double m2, double h)
{
	if (h < 0) 
		h += 1.0;
	if (h > 1) 
		h -= 1.0;
	if (6.0*h < 1)
		return (m1+(m2-m1)*h*6.0);
	if (2.0*h < 1)
		return m2;
	if (3.0*h < 2.0)
		return (m1+(m2-m1)*((2.0/3.0)-h)*6.0);
	
	return m1;
}

COLORREF CImage::HLStoRGB(const double &H, const double &L, const double &S)
{
	double r,g,b;
	double m1, m2;

	if (S==0) 
	{
		r=g=b=L;
	} 
	else 
	{
		if (L <= 0.5)
			m2 = L*(1.0+S);
		else
		m2 = L+S-L*S;
		m1 = 2.0*L-m2;
		r = HuetoRGB(m1, m2, H+1.0/3.0);
		g = HuetoRGB(m1, m2, H);
		b =	HuetoRGB(m1, m2, H-1.0/3.0);
	}

	return RGB((BYTE)(r*255),(BYTE)(g*255),(BYTE)(b*255));
}

BOOL CImage::Negate()
{
	if (m_hDib == NULL)
		return FALSE;

	BOOL bSuccess = TRUE;
	for (int y = 0; bSuccess && y < Height(); y++)
	{
		for (int x = 0; bSuccess && x < Width(); x++)
		{
			COLORREF c;
			if (GetPixel(x, y, c))
				bSuccess = SetPixel(x, y, RGB(255-GetRValue(c),255-GetGValue(c),255-GetBValue(c)));
			else
				bSuccess = FALSE;
		}
	}

	return bSuccess;
}

BOOL CImage::Mirror()
{
	if (m_hDib == NULL)
		return FALSE;

	BOOL bSuccess = TRUE;

	int nMiddle = Width()/2;
  
	for (int x=nMiddle; bSuccess && x>0; x--)
	{
		for (int y=0; bSuccess && y<Height(); y++)
		{
			COLORREF c1;
			int x1 = nMiddle + (nMiddle - x);
			bSuccess = GetPixel(x1, y, c1);
        
			COLORREF c2;
			bSuccess = bSuccess && GetPixel(x, y, c2);
			bSuccess = bSuccess && SetPixel(x, y, c1);
			bSuccess = bSuccess && SetPixel(x1, y, c2);
		}
	}

	return bSuccess;
}

BOOL CImage::Draw(HDC hDC, POINT ptOriginDst, SIZE sizeDst, RECT rcSrc)
{
	if(m_pBMI == NULL) return FALSE;
//	if(m_hPalette != NULL) {
//		::SelectPalette(pDC->GetSafeHdc(), m_hPalette, TRUE);
//	}
	int nXSrc, nYSrc, nWidthSrc, nHeightSrc;
	if (rcSrc.left == 0 && rcSrc.right == 0 && rcSrc.top == 0 && rcSrc.bottom == 0)
	{
		nXSrc = nYSrc = 0;
		nWidthSrc = Width();
		nHeightSrc = Height();
	}
	else
	{
		nXSrc = rcSrc.left;
		nYSrc = rcSrc.top;
		nWidthSrc = abs(rcSrc.right - rcSrc.left);
		nHeightSrc = abs(rcSrc.bottom - rcSrc.top);		
	}

	if(!(nXSrc >= 0 && nWidthSrc + nXSrc <= Width()))printf("nXSrc >= 0 && nWidthSrc + nXSrc <= Width() wrong/n");
	if(!(nYSrc >= 0 && nHeightSrc + nYSrc  <= Height()))printf("nYSrc >= 0 && nHeightSrc + nYSrc  <= Height()/n");

	::SetStretchBltMode(hDC,COLORONCOLOR);
	::StretchDIBits(hDC, 
		            ptOriginDst.x, 
					ptOriginDst.y, 
					sizeDst.cx, 
					sizeDst.cy,
		            nXSrc, 
					nYSrc, 
					nWidthSrc, 
					nHeightSrc,
		            m_pPixel, 
					m_pBMI, 
					DIB_RGB_COLORS, 
					SRCCOPY);
	return TRUE;
}

// bad function!!
/*
BOOL CImage::SaveDCAsFile(HDC hDC, LPCTSTR lpszFilePath)
{
	FILE* f = fopen(lpszFilePath, "wb");
	CMDIFrameWnd* pFrame = (CMDIFrameWnd*)AfxGetMainWnd();
	CImageChildFrame* pChild = (CImageChildFrame*)pFrame->GetActiveFrame();

	Lock();

	BITMAPFILEHEADER BmpFileHdr;
	BITMAPINFOHEADER* pBmpInfoHdr = (LPBITMAPINFOHEADER) m_pBMI;
	LONG lColors, lColorTableSize, lBitsSize, lBISize;

	// Determine size.

	lColors = GetColorEntries(m_pBMI);
	lColorTableSize = lColors * sizeof(RGBQUAD);
	lBISize = sizeof(BITMAPINFOHEADER) + lColorTableSize;
	DWORD dwWidthBytes = (pBmpInfoHdr->biWidth*pBmpInfoHdr->biBitCount + 31) / 32 * 4;
	lBitsSize = dwWidthBytes * pBmpInfoHdr->biHeight;

	// Write the file header.

	BmpFileHdr.bfType = 0x4D42;		// 'BM'
	BmpFileHdr.bfSize = sizeof(BITMAPFILEHEADER) + lBISize + lBitsSize;
	BmpFileHdr.bfReserved1 = 0;
	BmpFileHdr.bfReserved2 = 0;
	BmpFileHdr.bfOffBits = sizeof(BITMAPFILEHEADER) + lBISize;

	fwrite(&BmpFileHdr, sizeof(BITMAPFILEHEADER), 1, f);

	// Write the BITMAPINFOHEADER and color table.

	fwrite(m_pBMI, lBISize, 1, f);

	// Write bits data.
	COLORREF value;
	BYTE* pPixel = new BYTE[lBitsSize];
	int step = Width()*Height()/100;

	for(int i = 0; i < Width(); i ++)
		for( int j = 0; j < Height(); j ++)
		{
			BYTE* pData = pPixel + (Height()-1-j)*m_lWidthBytes + i*3;  //THIS NEEDS TO BE TESTED

			value = ::GetPixel(hDC,i,j);
			pData[0] = GetBValue(value);
			pData[1] = GetGValue(value);
			pData[2] = GetRValue(value);

			if(((i-1)*Height()+j)%step == 0)
				pChild->m_wndMiniStatusBar.m_progress.StepIt();
		}
	pChild->m_wndMiniStatusBar.m_progress.SetPos(0);
		
	fwrite(pPixel, lBitsSize, 1, f);

	Unlock();

	fclose(f);

	return TRUE;
}
*/
BOOL CImage::Save(LPCTSTR lpszFilePath, LPBYTE pData)
{
	FILE* f = fopen(lpszFilePath, "wb+");

	Lock();

	BITMAPFILEHEADER BmpFileHdr;
	BITMAPINFOHEADER* pBmpInfoHdr = (LPBITMAPINFOHEADER) m_pBMI;
	LONG lColors, lColorTableSize, lBitsSize, lBISize;

	// Determine size.

	lColors = GetColorEntries(m_pBMI);
	lColorTableSize = lColors * sizeof(RGBQUAD);
	lBISize = sizeof(BITMAPINFOHEADER) + lColorTableSize;
	DWORD dwWidthBytes = (pBmpInfoHdr->biWidth*pBmpInfoHdr->biBitCount + 31) / 32 * 4;
	lBitsSize = dwWidthBytes * pBmpInfoHdr->biHeight;

	// Write the file header.

	BmpFileHdr.bfType = 0x4D42;		// 'BM'
	BmpFileHdr.bfSize = sizeof(BITMAPFILEHEADER) + lBISize + lBitsSize;
	BmpFileHdr.bfReserved1 = 0;
	BmpFileHdr.bfReserved2 = 0;
	BmpFileHdr.bfOffBits = sizeof(BITMAPFILEHEADER) + lBISize;

	fwrite(&BmpFileHdr, 1, sizeof(BITMAPFILEHEADER), f);

	// Write the BITMAPINFOHEADER and color table.

	fwrite(m_pBMI, 1, lBISize, f);

	// Write bits data.

	fwrite(pData, 1, lBitsSize, f);

	Unlock();

	fclose(f);

	return TRUE;
}

void CImage::Convert2Gray(CImage *pDestImg)
{
	assert(m_pBMI->bmiHeader.biBitCount == 24);

	int width = Width();
	int height = Height();
	BYTE byte;
	
	// Set all members.
	for(int y = 0; y < height; y++)
		for(int x = 0; x < width; x ++)
		{
			RGBTRIPLE rgbPixel = *((RGBTRIPLE*)(m_pPixel + y * m_lWidthBytes + 3 * x));
            // gray = r *0.114f + g*0.588f + b*0.298
			byte = (BYTE)(((int)rgbPixel.rgbtRed* 3735 + 
								 (int)rgbPixel.rgbtGreen * 19267 + 
								 (int)rgbPixel.rgbtBlue * 9765) 
								 /32767);
			*(pDestImg->m_pPixel + y * m_lWidthBytes + 3 * x) = byte;
			*(pDestImg->m_pPixel + y * m_lWidthBytes + 3 * x + 1) = byte;
			*(pDestImg->m_pPixel + y * m_lWidthBytes + 3 * x + 2) = byte;
		}
}
