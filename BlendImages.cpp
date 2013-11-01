///////////////////////////////////////////////////////////////////////////
//
// NAME
//  BlendImages.cpp -- blend together a set of overlapping images
//
// DESCRIPTION
//  This routine takes a collection of images aligned more or less horizontally
//  and stitches together a mosaic.
//
//  The images can be blended together any way you like, but I would recommend
//  using a soft halfway blend of the kind Steve presented in the first lecture.
//
//  Once you have blended the images together, you should crop the resulting
//  mosaic at the halfway points of the first and last image.  You should also
//  take out any accumulated vertical drift using an affine warp.
//  Lucas-Kanade Taylor series expansion of the registration error.
//
// SEE ALSO
//  BlendImages.h       longer description of parameters
//
// Copyright ?Richard Szeliski, 2001.  See Copyright.h for more details
// (modified for CSE455 Winter 2003)
//
///////////////////////////////////////////////////////////////////////////

#include "ImageLib/ImageLib.h"
#include "BlendImages.h"
#include <float.h>
#include <math.h>
#include <iostream>
using namespace std;

#define MAX(x,y) (((x) < (y)) ? (y) : (x))
#define MIN(x,y) (((x) < (y)) ? (x) : (y))

/* Return the closest integer to x, rounding up */
static int iround(double x) {
    if (x < 0.0) {
        return (int) (x - 0.5);
    } else {
        return (int) (x + 0.5);
    }
}

void ImageBoundingBox(CImage &image, CTransform3x3 &M, 
    int &min_x, int &min_y, int &max_x, int &max_y)
{
    // This is a useful helper function that you might choose to implement
    // takes an image, and a transform, and computes the bounding box of the
    // transformed image.
	CShape sh=image.Shape();
    int width=sh.width;
    int height=sh.height;
	CVector3 p1,p2,p3,p4;
	p1[0]=0;p1[1]=0;p1[2]=1;
	p2[0]=0;p2[1]=height;p2[2]=1;
	p3[0]=width;p3[1]=0;p3[2]=1;
	p4[0]=width;p4[1]=height;p4[2]=1;
	p1=M*p1;p2=M*p2;p3=M*p3;p4=M*p4;
	min_x=iround(MIN(p1[0]/p1[2],p2[0]/p2[2]));
	max_x=iround(MAX(p3[0]/p3[2],p4[0]/p4[2]));
	min_y=iround(MIN(p1[1]/p1[2],p3[1]/p3[2]));
	max_y=iround(MAX(p2[1]/p2[2],p4[1]/p4[2]));
}


/******************* TO DO *********************
* AccumulateBlend:
*	INPUT:
*		img: a new image to be added to acc
*		acc: portion of the accumulated image where img is to be added
*       M: the transformation mapping the input image 'img' into the output panorama 'acc'
*		blendWidth: width of the blending function (horizontal hat function;
*	    try other blending functions for extra credit)
*	OUTPUT:
*		add a weighted copy of img to the subimage specified in acc
*		the first 3 band of acc records the weighted sum of pixel colors
*		the fourth band of acc records the sum of weight
*/
static void AccumulateBlend(CByteImage& img, CFloatImage& acc, CTransform3x3 M, float blendWidth)
{
    // BEGIN TODO
    // Fill in this routine
	// get shape of acc and img
	CShape sh=img.Shape();
    int width=sh.width;
    int height=sh.height;
	CShape shacc=acc.Shape();
    int widthacc=shacc.width;
    int heightacc=shacc.height;
	
	// get the bounding box of img in acc
	int min_x,min_y,max_x,max_y;
	ImageBoundingBox(img,M,min_x,min_y,max_x,max_y);
	int middle_x=(max_x+min_x)/2;
	
	// --> Trying another way to do exposure compensation
	/*for (int ii = 0; ii < widthacc; ii++)
		for (int jj = 0; jj < heightacc; jj++)
		{
			CVector3 p;
			p[0] = ii; p[1] = jj; p[2] = 1;
			p = M.Inverse() * p;
			int newx, newy;
			newx = iround(p[0] / p[2]);
			newy = iround(p[1] / p[2]);
			if (acc.Pixel(ii,jj,0) != 0 && acc.Pixel(ii,jj,1) != 0 && acc.Pixel(ii,jj,2) != 0 && img.Pixel(;
		}*/
	// <-- Trying another way to do exposure compensation

	// add every pixel in img to acc, feather the region withing blendwidth to the bounding box,
	// pure black pixels (caused by warping) are not added
	CVector3 p;
	//int newx,newy;
	double weight;
	/*for (int ii=0;ii<width;ii++)
	{
        for (int jj=0;jj<height;jj++)
		{
		    p[0]=ii;p[1]=jj;p[2]=1;
			p=M*p;
			//newx=iround(p[0]/p[2]);newy=iround(p[1]/p[2]);
			if ((newx>=0)&&(newx<widthacc)&&(newy>=0)&&(newy<heightacc))
			{
				weight=1;
				if ((newx>min_x)&&(newx<(min_x+blendWidth)))
					weight=(newx-min_x)/blendWidth;
				if ((newx<max_x)&&(newx>(max_x-blendWidth)))
					weight=(max_x-newx)/blendWidth;
				if ((img.Pixel(ii,jj,0)==0)&&(img.Pixel(ii,jj,0)==0)&&(img.Pixel(ii,jj,0)==0))
					weight=0;
				acc.Pixel(newx,newy,0)+=(img.Pixel(ii,jj,0)*weight);
				acc.Pixel(newx,newy,1)+=(img.Pixel(ii,jj,1)*weight);
				acc.Pixel(newx,newy,2)+=(img.Pixel(ii,jj,2)*weight);
				acc.Pixel(newx,newy,3)+=weight;
			}
		}
	}*/
	
	double newx, newy;
	// --> Trying E C
	float diff_r = 0.0;
	float diff_g = 0.0;
	float diff_b = 0.0;
	int cnt = 0;
	for (int ii = 0; ii < widthacc; ii++)
	{
		for (int jj = 0; jj < heightacc; jj++)
		{
			p[0] = ii; p[1] = jj; p[2] = 1;
			p = M.Inverse() * p;
			newx = p[0] / p[2];
			newy = p[1] / p[2];
			bool flag = true;
			if (newx>=0 && newx<width-1 && newy>=0 && newy<height-1)
			{
				int newxx = (int)newx;
				int newyy = (int)newy;
				if (acc.Pixel(ii,jj,0)==0
					&& acc.Pixel(ii,jj,1)==0
					&& acc.Pixel(ii,jj,2)==0)
					flag = false;
				if (img.Pixel(newxx,newyy,0)==0
					&& img.Pixel(newxx,newyy,1)==0
					&& img.Pixel(newxx,newyy,2)==0)
					flag = false;
				if (img.Pixel(newxx+1,newyy,0)==0
					&& img.Pixel(newxx+1,newyy,1)==0
					&& img.Pixel(newxx+1,newyy,2)==0)
					flag = false;
				if (img.Pixel(newxx,newyy+1,0)==0
					&& img.Pixel(newxx,newyy+1,1)==0
					&& img.Pixel(newxx,newyy+1,2)==0)
					flag = false;
				if (img.Pixel(newxx+1,newyy+1,0)==0
					&& img.Pixel(newxx+1,newyy+1,1)==0
					&& img.Pixel(newxx+1,newyy+1,2)==0)
					flag = false;
				if (flag)
				{
					cnt++;
					diff_r += acc.Pixel(ii,jj,0) - img.PixelLerp(newx,newy,0);
					diff_g += acc.Pixel(ii,jj,1) - img.PixelLerp(newx,newy,1);
					diff_b += acc.Pixel(ii,jj,2) - img.PixelLerp(newx,newy,2);
				}
				if (cnt > 10 && cnt <= 20)
					cout << "cnt" << cnt << " " << diff_r << " " << diff_g << " " << diff_b << endl;
			}
		}
	}
	if (cnt != 0)
	{
		diff_r /= (float) cnt;
		diff_g /= (float) cnt;
		diff_b /= (float) cnt;
	}
	cout << "Compensation: " << diff_r << " " << diff_g << " " << diff_b << endl;
	// <-- Trying E C
	for (int ii = 0; ii < widthacc; ii++)
		for (int jj = 0; jj < heightacc; jj++)
		{
			p[0] = ii; p[1] = jj; p[2] = 1;
			p = M.Inverse() * p;
			newx = p[0] / p[2];
			newy = p[1] / p[2];
			if ((newx>=0)&&(newx<width)&&(newy>=0)&&(newy<height))
			{
				weight=1;
				if ((ii>min_x)&&(ii<(min_x+blendWidth)))
					weight=(ii-min_x)/blendWidth;
				if ((ii<max_x)&&(ii>(max_x-blendWidth)))
					weight=(max_x-ii)/blendWidth;
				if ((img.Pixel(newx,newy,0)==0)&&(img.Pixel(newx,newy,1)==0)&&(img.Pixel(newx,newy,2)==0))
					weight=0;
				acc.Pixel(ii,jj,0)+=((img.PixelLerp(newx,newy,0)+diff_r)*weight);
				acc.Pixel(ii,jj,1)+=((img.PixelLerp(newx,newy,1)+diff_g)*weight);
				acc.Pixel(ii,jj,2)+=((img.PixelLerp(newx,newy,2)+diff_b)*weight);
				acc.Pixel(ii,jj,3)+=weight;
			}
		}
	
	printf("AccumulateBlend\n"); 

    // END TODO
}



/******************* TO DO 5 *********************
* NormalizeBlend:
*	INPUT:
*		acc: input image whose alpha channel (4th channel) contains
*		     normalizing weight values
*		img: where output image will be stored
*	OUTPUT:
*		normalize r,g,b values (first 3 channels) of acc and store it into img
*/
static void NormalizeBlend(CFloatImage& acc, CByteImage& img)
{
    // BEGIN TODO
    // fill in this routine..
	// divide the total weight for every pixel
	CShape shacc=acc.Shape();
    int widthacc=shacc.width;
    int heightacc=shacc.height;
	for (int ii=0;ii<widthacc;ii++)
	{
        for (int jj=0;jj<heightacc;jj++)
		{
			if (acc.Pixel(ii,jj,3)>0)
			{
				img.Pixel(ii,jj,0)=(int)(acc.Pixel(ii,jj,0)/acc.Pixel(ii,jj,3));
				img.Pixel(ii,jj,1)=(int)(acc.Pixel(ii,jj,1)/acc.Pixel(ii,jj,3));
				img.Pixel(ii,jj,2)=(int)(acc.Pixel(ii,jj,2)/acc.Pixel(ii,jj,3));
			}
			else
			{
				img.Pixel(ii,jj,0)=0;
				img.Pixel(ii,jj,1)=0;
				img.Pixel(ii,jj,2)=0;
			}
		}
	}

    // END TODO
}



/******************* TO DO 5 *********************
* BlendImages:
*	INPUT:
*		ipv: list of input images and their relative positions in the mosaic
*		blendWidth: width of the blending function
*	OUTPUT:
*		create & return final mosaic by blending all images
*		and correcting for any vertical drift
*/
CByteImage BlendImages(CImagePositionV& ipv, float blendWidth)
{
    // Assume all the images are of the same shape (for now)
    CByteImage& img0 = ipv[0].img;
    CShape sh        = img0.Shape();
    int width        = sh.width;
    int height       = sh.height;
    int nBands       = sh.nBands;
    // int dim[2]       = {width, height};

    int n = ipv.size();
    if (n == 0) return CByteImage(0,0,1);

    bool is360 = false;

    // Hack to detect if this is a 360 panorama
    if (ipv[0].imgName == ipv[n-1].imgName)
        is360 = true;

    // Compute the bounding box for the mosaic
    float min_x = FLT_MAX, min_y = FLT_MAX;
    float max_x = 0, max_y = 0;
    int i;
	int tmpmin_x,tmpmin_y,tmpmax_x,tmpmax_y;
    for (i = 0; i < n; i++)
    {
        CTransform3x3 &T = ipv[i].position;

        // BEGIN TODO
        // add some code here to update min_x, ..., max_y
		ImageBoundingBox(img0,T,tmpmin_x,tmpmin_y,tmpmax_x,tmpmax_y);
		min_x=MIN(tmpmin_x,min_x);
		min_y=MIN(tmpmin_y,min_y);
		max_x=MAX(tmpmax_x,max_x);
		max_y=MAX(tmpmax_y,max_y);
		// END TODO
    }

    // Create a floating point accumulation image
    CShape mShape((int)(ceil(max_x) - floor(min_x)), (int)(ceil(max_y) - floor(min_y)), nBands + 1);
    CFloatImage accumulator(mShape);
    accumulator.ClearPixels();

    double x_init, x_final;
    double y_init, y_final;

	//--> Exposure Compensation
	/*float* alphaR = new float[n];
	float* alphaG = new float[n];
	float* alphaB = new float[n];
	alphaR[0] = 1; alphaG[0] = 1; alphaB[0] = 1;
	int minX_prev, minY_prev, maxX_prev, maxY_prev;
	int minX_curr, minY_curr, maxX_curr, maxY_curr;
	// Computing correction coefficient alpha
	for (i = 1; i < n; i++)
	{
		float sum_prev[3]; // rgb
		float sum_curr[3]; // rgb
		CTransform3x3 M_prev = CTransform3x3::Translation(-min_x,-min_y)*ipv[i-1].position;
		CTransform3x3 M_curr = CTransform3x3::Translation(-min_x,-min_y)*ipv[i].position;
		ImageBoundingBox(ipv[i-1].img,M_prev,minX_prev,minY_prev,maxX_prev,maxY_prev);
		ImageBoundingBox(ipv[i].img,M_curr,minX_curr,minY_curr,maxX_curr,maxY_curr);
		int OMinX = minX_curr > minX_prev ? MIN(maxX_prev,minX_curr) : MIN(minX_prev,maxX_curr);
		int OMaxX = minX_curr > minX_prev ? MAX(maxX_prev,minX_curr) : MAX(minX_prev,maxX_curr);
		int OMinY = MAX(minY_prev,minY_curr);
		int OMaxY = MIN(maxY_prev,maxY_curr);
		CVector3 p,p_prev,p_curr;
		for (int m = OMinX; m <= OMaxX; m++)
			for (int n = OMinY; n <= OMaxY; n++)
			{
				p[0] = m; p[1] = n; p[2] = 1;
				p_prev = M_prev.Inverse() * p;
				p_curr = M_curr.Inverse() * p;
				for (int ch = 0; ch < 3; ch++)
				{
					sum_prev[ch] += (float)ipv[i-1].img.Pixel(p_prev[0]/p_prev[2],p_prev[1]/p_prev[2],ch) / 255.0;
					sum_curr[ch] += (float)ipv[i].img.Pixel(p_curr[0]/p_curr[2],p_curr[1]/p_curr[2],ch) / 255.0;
				}
			}
		alphaR[i] = sum_prev[0] / sum_curr[0];
		alphaG[i] = sum_prev[1] / sum_curr[1];
		alphaB[i] = sum_prev[2] / sum_curr[2];
	}
	// Computing global compensation coefficient g
	float gR, gG, gB;
	float sum_alphaR = 0.0, sum_alphaR2 = 0.0;
	float sum_alphaG = 0.0, sum_alphaG2 = 0.0;
	float sum_alphaB = 0.0, sum_alphaB2 = 0.0;
	for (i = 0; i < n; i++)
	{
		sum_alphaR += alphaR[i];
		sum_alphaR2 += alphaR[i] * alphaR[i];
		sum_alphaG += alphaG[i];
		sum_alphaG2 += alphaG[i] * alphaG[i];
		sum_alphaB += alphaB[i];
		sum_alphaB2 += alphaB[i] * alphaB[i];
	}
	gR = sum_alphaR / sum_alphaR2;
	gG = sum_alphaG / sum_alphaG2;
	gB = sum_alphaB / sum_alphaB2;*/
	// Compensation
	/*for (i = 0; i < n; i++)
	{
		CByteImage& imgi = ipv[i].img;
		int widthi = imgi.Shape().width;
		int heighti = imgi.Shape().height;
		float gamma1 = 1.0/2.2;
		for (int m = 0; m < widthi; m++)
			for (int n = 0; n < heighti; n++)
			{
				imgi.Pixel(m,n,0) *= pow(gR*alphaR[i],gamma1);
				cout << (float)imgi.Pixel(m,n,0) << endl;
				imgi.Pixel(m,n,1) *= pow(gG*alphaG[i],gamma1);
				imgi.Pixel(m,n,2) *= pow(gB*alphaB[i],gamma1);
			}
	}*/
	//<-- Exposure Compensation

    // Add in all of the images
    for (i = 0; i < n; i++) 
	{
        // Compute the sub-image involved
        CTransform3x3 &M = ipv[i].position;
        CTransform3x3 M_t = CTransform3x3::Translation(-min_x, -min_y) * M;
        CByteImage& img = ipv[i].img;
		/*float gamma1 = 1.0/2.2;
		for (int m = 0; m < img.Shape().width; m++)
			for (int n = 0; n < img.Shape().height; n++)
			{
				img.Pixel(m,n,0) = floor((float)img.Pixel(m,n,0)*gR*alphaR[i]);//pow(gR*alphaR[i],gamma1));
				//cout << (float)imgi.Pixel(m,n,0) << endl;
				img.Pixel(m,n,1) = floor((float)img.Pixel(m,n,1)*gG*alphaG[i]);//*pow(gG*alphaG[i],gamma1));
				img.Pixel(m,n,2) = floor((float)img.Pixel(m,n,2)*gB*alphaB[i]);//*pow(gB*alphaB[i],gamma1));
			}
		if (i == 0)
			WriteFile(img, "tmp_img0.tga");
		else if (i == 1)
			WriteFile(img, "tmp_img1.tga");
		else if (i == 2)
			WriteFile(img, "tmp_img2.tga");
		else if (i == 3)
			WriteFile(ipv[0].img, "ref_img0.tga");*/
		// --> Added for exposure compensation
		/*int minX1, minY1, maxX1, maxY1;
		int minX2, minY2, maxX2, maxY2;
		if (i == 0)
			ImageBoundingBox(img, M_t, minX1, minY1, maxX1, maxY1);
		else if (i == 1)
			ImageBoundingBox(img, M_t, minX2, minY2, maxX2, maxY2);
		else
		{
			minX1 = minX2;
			minY1 = minY2;
			maxX1 = maxX2;
			maxY1 = maxY2;
			ImageBoundingBox(img, M_t, minX2, minY2, maxX2, maxY2);
		}
		if (i > 0 && i < n)
		{
			CByteImage& imgprev = ipv[i-1].img;
			CTransform3x3 Mprev = CTransform3x3::Translation(-min_x, -min_y) * ipv[i-1].position;
			int OverlapMinX = minX2 > minX1 ? MIN(maxX1,minX2) : MIN(minX1,maxX2);
			int OverlapMaxX = minX2 > minX1 ? MAX(maxX1,minX2) : MAX(minX1,maxX2);;
			int OverlapMinY = MAX(minY1,minY2);
			int OverlapMaxY = MIN(maxY1,maxY2);
			float diff_r = 0.0, diff_g = 0.0, diff_b = 0.0;
			int cnt = 0;
			for (int m = OverlapMinX; m <= OverlapMaxX; m++)
			{
				for (int n = OverlapMinY; n <= OverlapMaxY; n++)
				{
					CVector3 p,q;
					p[0] = m;
					p[1] = n;
					p[2] = 1;
					q = p;
					p = M_t.Inverse()*p;
					q = Mprev.Inverse()*q;
					int mm = p[0]/p[2];
					int nn = p[1]/p[2];
					int mmm = q[0]/q[2];
					int nnn = q[1]/q[2];
					if (imgprev.Pixel(mmm,nnn,0)!=0 && imgprev.Pixel(mmm,nnn,1)!=0 &&
						imgprev.Pixel(mmm,nnn,2)!=0 && img.Pixel(mm,nn,0)!=0 &&
						img.Pixel(mm,nn,1)!=0 && img.Pixel(mm,nn,2)!=0)
					{
						cnt++;
						diff_r += imgprev.Pixel(mmm,nn,0) - img.Pixel(mm,nn,0);
						diff_g += imgprev.Pixel(mmm,nn,1) - img.Pixel(mm,nn,1);
						diff_b += imgprev.Pixel(mmm,nn,2) - img.Pixel(mm,nn,2);
						if (cnt <= 10)
							cout << "cnt: " << cnt << " " << diff_r << " " << diff_g << " " << diff_b << endl;
						if (cnt == 10)
						{
							cout << OverlapMinX << " " << OverlapMaxX << endl;
							cout << OverlapMinY << " " << OverlapMaxY << endl;
						}
					}
				}
			}
			cout << "diffsum_r: " << diff_r << endl;
			cout << "diffsum_g: " << diff_g << endl;
			cout << "diffsum_b: " << diff_b << endl;
			diff_r /= (float) cnt;
			diff_g /= (float) cnt;
			diff_b /= (float) cnt;
			cout << "diff_r: " << diff_r << endl;
			cout << "diff_g: " << diff_g << endl;
			cout << "diff_b: " << diff_b << endl;
			for (int m = 0; m <= width; m++)
				for (int n = 0; n <= height; n++)
				{
					if (img.Pixel(m,n,0)!=0 && img.Pixel(m,n,1)!=0 && img.Pixel(m,n,2)!=0)
					{
						img.Pixel(m,n,0) += diff_r;
						img.Pixel(m,n,1) += diff_g;
						img.Pixel(m,n,2) += diff_b;
						if (img.Pixel(m,n,0) > 255) img.Pixel(m,n,0) = 255;
						if (img.Pixel(m,n,1) > 255) img.Pixel(m,n,1) = 255;
						if (img.Pixel(m,n,2) > 255) img.Pixel(m,n,2) = 255;
						if (img.Pixel(m,n,0) < 0) img.Pixel(m,n,0) = 0;
						if (img.Pixel(m,n,1) < 0) img.Pixel(m,n,1) = 0;
						if (img.Pixel(m,n,2) < 0) img.Pixel(m,n,2) = 0;
					}
				}
		}*/
		// <-- Added for exposure compensation

        // Perform the accumulation
        AccumulateBlend(img, accumulator, M_t, blendWidth);
		
        if (i == 0) 
		{
            CVector3 p;
            p[0] = 0.5 * width;
            p[1] = 0.0;
            p[2] = 1.0;

            p = M_t * p;
            x_init = p[0];
            y_init = p[1];
        } else if (i == n - 1) 
		{
            CVector3 p;
            p[0] = 0.5 * width;
            p[1] = 0.0;
            p[2] = 1.0;

            p = M_t * p;
            x_final = p[0];
            y_final = p[1];
        }
    }

    // Normalize the results
    mShape = CShape((int)(ceil(max_x) - floor(min_x)), (int)(ceil(max_y) - floor(min_y)), nBands);

    CByteImage compImage(mShape);
    NormalizeBlend(accumulator, compImage);
    bool debug_comp = false;
    if (debug_comp)
        WriteFile(compImage, "tmp_comp.tga");

    // Allocate the final image shape
    int outputWidth = 0;
    if (is360) 
	{
        outputWidth = mShape.width - width;
    } else 
	{
        outputWidth = mShape.width;
    }

    CShape cShape(outputWidth, mShape.height, nBands);

    CByteImage croppedImage(cShape);

    // Compute the affine transformation
    CTransform3x3 A = CTransform3x3(); // identify transform to initialize

    // BEGIN TODO
    // fill in appropriate entries in A to trim the left edge and
    // to take out the vertical drift if this is a 360 panorama
    // (i.e. is360 is true)
	double k;
	if (is360)
	{
		if (x_init>x_final)
		{
			int tmp=x_init;
			x_init=x_final;
			x_final=tmp;
		}
		CTransform3x3 AA = CTransform3x3();;
		k=-(y_final-y_init)/(x_final-x_init);
		AA[0][0]=1;AA[0][1]=0;AA[0][2]=0;
		AA[1][0]=k;AA[1][1]=1;AA[1][2]=0;
		AA[2][0]=0;AA[2][1]=0;AA[2][2]=1;
		A = CTransform3x3::Translation(x_init, -y_init) * AA;
	}

    // END TODO

    // Warp and crop the composite
    WarpGlobal(compImage, croppedImage, A, eWarpInterpLinear);

    return croppedImage;
}

