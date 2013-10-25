///////////////////////////////////////////////////////////////////////////
//
// NAME
//  FeatureAlign.cpp -- image registration using feature matching
//
// SEE ALSO
//  FeatureAlign.h      longer description
//
// Based on code by Richard Szeliski, 2001.
// (modified for CSE576 Spring 2005, and for CS4670, Fall 2012-2013)
//
///////////////////////////////////////////////////////////////////////////

#include "ImageLib/ImageLib.h"
#include "FeatureAlign.h"
#include "SVD.h"

#include <math.h>
#include <time.h>
#include <iostream>

CTransform3x3 ComputeHomography(const FeatureSet &f1, const FeatureSet &f2,
                                const vector<FeatureMatch> &matches)
{
    int numMatches = (int) matches.size();

    // first, we will compute the A matrix in the homogeneous linear equations Ah = 0
    int numRows = 2 * numMatches; // number of rows of A
    const int numCols = 9;        // number of columns of A

    // this allocates space for the A matrix
    AMatrixType A = AMatrixType::Zero(numRows, numCols);

    for (int i = 0; i < numMatches; i++) {
        const FeatureMatch &m = matches[i];
        const Feature &a = f1[m.id1];
        const Feature &b = f2[m.id2];

        // BEGIN TODO
        // fill in the matrix A in this loop.
        // To access an element of A, use parentheses, e.g. A(0,0)

		// Note that A is a (2n x 9) matrix
		// We can easily fill in A according to the lecture
		// Row 2i
		A(2*i, 0) = a.x;
		A(2*i, 1) = a.y;
		A(2*i, 2) = 1;
		A(2*i, 3) = 0;
		A(2*i, 4) = 0;
		A(2*i, 5) = 0;
		A(2*i, 6) = -b.x * a.x;
		A(2*i, 7) = -b.x * a.y;
		A(2*i, 8) = -b.x;

		// Row 2i+1
		A(2*i+1, 0) = 0;
		A(2*i+1, 1) = 0;
		A(2*i+1, 2) = 0;
		A(2*i+1, 3) = a.x;
		A(2*i+1, 4) = a.y;
		A(2*i+1, 5) = 1;
		A(2*i+1, 6) = -b.y * a.x;
		A(2*i+1, 7) = -b.y * a.y;
		A(2*i+1, 8) = -b.y;

        // END TODO
    }

    // Compute the SVD of the A matrix and get out the matrix V^T and the vector of singular values
    AMatrixType Vt;
    VectorXd sv;
    SVD(A, Vt, sv);

    CTransform3x3 H;
    // BEGIN TODO
    // fill the homography H with the appropriate elements of the SVD
    // To extract, for instance, the V matrix, use svd.matrixV()

	// Find the smallest sv first
	// minSv: smallest sv value
	// minSvIdx: smallest sv index
	double minSv = sv[0];
	int minSvIdx = 0;
	for (int i = 1; i < (int)sv.size(); i++) {
		if (sv[i] < minSv && sv[i] > 0) {
			minSv = sv[i];
			minSvIdx = i;
		}
	}
	// Rank(A) = 8 => H = a v_n  (a is a constant)
	// v_n corresponds to the smallest sv
	H[0][0] = Vt(minSvIdx,0) / Vt(minSvIdx,8);
	H[0][1] = Vt(minSvIdx,1) / Vt(minSvIdx,8);
	H[0][2] = Vt(minSvIdx,2) / Vt(minSvIdx,8);
	H[1][0] = Vt(minSvIdx,3) / Vt(minSvIdx,8);
	H[1][1] = Vt(minSvIdx,4) / Vt(minSvIdx,8);
	H[1][2] = Vt(minSvIdx,5) / Vt(minSvIdx,8);
	H[2][0] = Vt(minSvIdx,6) / Vt(minSvIdx,8);
	H[2][1] = Vt(minSvIdx,7) / Vt(minSvIdx,8);
	H[2][2] = Vt(minSvIdx,8) / Vt(minSvIdx,8);

    // END TODO

    return H;
}


/******************* TO DO *********************
 * alignPair:
 *	INPUT:
 *		f1, f2: source feature sets
 *		matches: correspondences between f1 and f2
 *               Each match in 'matches' contains two feature ids of 
 *               matching features, id1 (in f1) and id2 (in f2).
 *		m: motion model
 *		nRANSAC: number of RANSAC iterations
 *		RANSACthresh: RANSAC distance threshold
 *		M: transformation matrix (output)
 *
 *	OUTPUT:
 *		repeat for nRANSAC iterations:
 *			choose a minimal set of feature matches
 *			estimate the transformation implied by these matches
 *			count the number of inliers
 *		for the transformation with the maximum number of inliers,
 *		compute the least squares motion estimate using the inliers,
 *		and store it in M
 */
int alignPair(const FeatureSet &f1, const FeatureSet &f2,
          const vector<FeatureMatch> &matches, MotionModel m, 
          int nRANSAC, double RANSACthresh, CTransform3x3& M)
{
    // BEGIN TODO
    // Write this entire method.  You need to handle two types of 
    // motion models, pure translations (m == eTranslation) and 
    // full homographies (m == eHomography).  However, you should
    // only have one outer loop to perform the RANSAC code, as 
    // the use of RANSAC is almost identical for both cases.
    //
    // Your homography handling code should call ComputeHomography.
    // This function should also call countInliers and, at the end,
    // leastSquaresFit.
	
	// Generate a random seed
	srand(unsigned int(time(NULL)));

	// Save maximum inliers
	// maxInlierCnt: maximum number of inliers
	// maxInliers: inlier IDs
	int maxInlierCnt = 0;
	vector<int> maxInliers;

	// Number of features
	int numF = matches.size();

	// Number of samples
	int numS;
	switch (m) {
		case eTranslate: {
			numS = 1;
			break;
		}
		case eHomography: {
			numS = 4;
			break;
		}
	}

	for (int i = 0; i < nRANSAC; i++) {
		// Random samples
		vector<FeatureMatch> randMatches;
		randMatches.clear();
		
		int cnt = 0;
		while (cnt < numS) {
			// Randomly pick a pair from the feature list that isn't in the list
			int rnd = rand() % numF;
			FeatureMatch curFeature = matches[rnd];
			// Check if this match is already selected
			bool selected = false;
			for (int j = 0; j < (int)randMatches.size(); j++) {
				if (randMatches[j].id1 == curFeature.id1 && randMatches[j].id2 == curFeature.id2) {
					selected = true;
				}
			}
			// If not, push it to the sample space
			if (!selected) {
				randMatches.push_back(curFeature);
				cnt++;
			}
		}
		
		// Get the transform
		CTransform3x3 H;
		// Translation
		if (numS == 1) {
			FeatureMatch match = randMatches[0];
			H[0][0] = 1;
			H[0][1] = 0;
			H[0][2] = f2[match.id2].x - f1[match.id1].x;
			H[1][0] = 0;
			H[1][1] = 1;
			H[1][2] = f2[match.id2].y - f1[match.id1].y;
			H[2][0] = 0;
			H[2][1] = 0;
			H[2][2] = 1;
		}
		// Homography
		else if (numS == 4) {
			H = ComputeHomography(f1, f2, randMatches);
		}

		// Count inliers matching this homography
		vector<int> curInliers;
		int numInliers = countInliers(f1, f2, matches, m, H, RANSACthresh, curInliers);

		// Update maximum inliers
		if (numInliers > maxInlierCnt) {
			maxInlierCnt = numInliers;
			maxInliers.clear();
			for (int k = 0; k < maxInlierCnt; k++) {
				maxInliers.push_back(curInliers[k]);
			}
		}
	}
	cout << "num_inliers: " << maxInlierCnt << " / " << matches.size() << endl;
	
	//////////////////////////DEBUG//////////////////////////
	/*if (maxInlierCnt < 10) {
		for (int ii = 0; ii < maxInlierCnt; ii++) {
			cout << ii << endl;
			cout << matches[maxInliers[ii]].id1 << " " << matches[maxInliers[ii]].id2 << endl;
			cout << f1[matches[maxInliers[ii]].id1].x << " " << f1[matches[maxInliers[ii]].id1].y << endl;
			cout << f2[matches[maxInliers[ii]].id2].x << " " << f2[matches[maxInliers[ii]].id2].y << endl;
		}
	}*/

	// Call leastSquaresFit
	leastSquaresFit(f1, f2, matches, m, maxInliers, M);
	
    // END TODO

    return 0;
}

/******************* TO DO *********************
 * countInliers:
 *	INPUT:
 *		f1, f2: source feature sets
 *		matches: correspondences between f1 and f2
 *		m: motion model
 *               Each match in 'matches' contains two feature ids of 
 *               matching features, id1 (in f1) and id2 (in f2).
 *		M: transformation matrix
 *		RANSACthresh: RANSAC distance threshold
 *		inliers: inlier feature IDs
 *	OUTPUT:
 *		transform the features in f1 by M
 *
 *		count the number of features in f1 for which the transformed
 *		feature is within Euclidean distance RANSACthresh of its match
 *		in f2
 *
 *		store these features IDs in inliers
 *
 */
int countInliers(const FeatureSet &f1, const FeatureSet &f2,
                 const vector<FeatureMatch> &matches, MotionModel m, 
                 CTransform3x3 M, double RANSACthresh, vector<int> &inliers)
{
    inliers.clear();

    for (unsigned int i = 0; i < matches.size(); i++) {
        // BEGIN TODO
        // determine if the ith matched feature f1[id1-1], when transformed by M,
        // is within RANSACthresh of its match in f2
        //
        // if so, append i to inliers

		// Get the features
		const FeatureMatch &m = matches[i];
		const Feature &a = f1[m.id1];
		const Feature &b = f2[m.id2];

		// Transform f1[m.id1] by M
		CVector3 p;
		p[0] = a.x;
		p[1] = a.y;
		p[2] = 1;
		CVector3 pt = M * p;
		
		// Transformed f1[m.id1]
		double xt = pt[0];
		double yt = pt[1];

		// Compute Euclidean distance from xt,yt to b.x,b.y
		// and check if within RANSACthresh
		double distance = sqrt((b.x-xt)*(b.x-xt) + (b.y-yt)*(b.y-yt));
		if (distance <= RANSACthresh) {
			// Append i to inliers
			inliers.push_back(i);
		}

        // END TODO
    }

    return (int) inliers.size();
}

/******************* TO DO *********************
 * leastSquaresFit:
 *	INPUT:
 *		f1, f2: source feature sets
 *		matches: correspondences between f1 and f2
 *		m: motion model
 *      inliers: inlier match indices (indexes into 'matches' array)
 *		M: transformation matrix (output)
 *	OUTPUT:
 *		compute the transformation from f1 to f2 using only the inliers
 *		and return it in M
 */
int leastSquaresFit(const FeatureSet &f1, const FeatureSet &f2,
            const vector<FeatureMatch> &matches, MotionModel m, 
            const vector<int> &inliers, CTransform3x3& M)
{
    // This function needs to handle two possible motion models, 
    // pure translations and full homographies.

    switch (m) {
        case eTranslate: {
            // for spherically warped images, the transformation is a 
            // translation and only has two degrees of freedom
            //
            // therefore, we simply compute the average translation vector
            // between the feature in f1 and its match in f2 for all inliers
            double u = 0;
            double v = 0;

            for (int i=0; i < (int) inliers.size(); i++) {
			    // BEGIN TODO
			    // use this loop to compute the average translation vector
			    // over all inliers
				FeatureMatch m = matches[inliers[i]];
				Feature a = f1[m.id1];
				Feature b = f2[m.id2];

				// Sum the differences
				u += b.x - a.x;
				v += b.y - a.y;

                // END TODO
            }

            u /= inliers.size();
            v /= inliers.size();

            M = CTransform3x3::Translation((float) u, (float) v);

            break;
        } 

        case eHomography: {
			M = CTransform3x3();

            // BEGIN TODO
		    // Compute a homography M using all inliers.
		    // This should call ComputeHomography.

			// Put the matches specified by inliers into a new FeatureMatch vector
			vector<FeatureMatch> inlierMatches;
			inlierMatches.clear();
			for (int i = 0; i < (int)inliers.size(); i++) {
				FeatureMatch m = matches[inliers[i]];
				inlierMatches.push_back(m);
			}

			// Compute the homography using all inliers
			M = ComputeHomography(f1, f2, inlierMatches);

            // END TODO
        
            break;
        }
    }    

    return 0;
}

