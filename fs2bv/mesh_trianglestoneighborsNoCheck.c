/*

using triangles of mesh to create neighbors list

FORMAT:       [nei, bn [, trl]] = meshtrianglestoneighbors(nc, tri)

Input fields:

      nc          1x1 number of vertices
      tri         Cx3 triangles (one-based)

Output fields:

      nei         Nx2 neighbors list (1-based)
      bn          1xN neighbors list with warnings
      tri         Nx1 triangle list (1-based)


% Version:  v0.7f
% Build:    9021300
% Date:     Feb-13 2009, 12:37 AM CET
% Author:   Jochen Weber, SCAN Unit, Columbia University, NYC, NY, USA
% URL/Info: http://wiki.brainvoyager.com/BVQXtools
%
% Modified from nesh_trianglestoneighbors.c. This function computes
% neighbors of the target nodes, ignoring any errors and warnings for
% bad nodes in which the neighbors are not assigned properly. This
% "ignoring" step is sometimes required to convert FreeSurfer-processed
% surface into BrainVoyager SRF files sicne the definitions of the "closed"
% surface seems to be different between those two software packages.
%
% Last Update: "2018-05-24 13:11:17 ban"

*/

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include <stdio.h>

#define MAXNRNEI 30
#define MAXNRNE2 60
#define MAXNRNEM 59

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    /* number of coordinates and counter (including double + triple) */
    double dnc, dn1, dn2, dn3;
    double *dnei, *snei;
    double dneib[MAXNRNEI], dneis[MAXNRNEI];
    int nc, cc, nt, nt2, sc;

    /* neighbors list pointer */
    int *nei, *neip;
    unsigned char *nnei, *nneic;
    unsigned char cinei, cnnei, tnnei;
    int n1, n2, n3;
    mxArray *cnei;
    mxArray *onnei;
    mxArray *onei;

    /* double array for bad neighbors */
    bool trackbn = false;
    double *bnd = NULL;
    signed long bnc = 0;
    unsigned long *bnl = NULL;

    /* number of dimensions and dim size */
    int nd;
    const int *di;
    int cdi[2] = {1, 2};
    int ndi[2] = {1, 1};
    
    int lcnt;

    /* variable output string */
    char vstr[256];

    /* check number of in/out arguments */
    if ((nrhs != 2) | (nlhs > 1))
            mexErrMsgTxt("Bad number of input/output arguments.");

    /* check argument types */
    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || (mxGetNumberOfElements(prhs[0]) != 1))
        mexErrMsgTxt("Input arguments must be of type double.");
    dnc = *((double *) mxGetPr(prhs[0]));

    /* good number */
    if (mxIsInf(dnc) || mxIsNaN(dnc) || dnc < 1 || dnc > 10000000 || dnc != ((double) ((int) dnc)))
        mexErrMsgTxt("Bad number of coordinates given.");
    nc = (int) dnc;

    /* check dims */
    nd = mxGetNumberOfDimensions(prhs[1]);
    di = mxGetDimensions(prhs[1]);
    if (nd != 2 || di[1] != 3)
        mexErrMsgTxt("Bad size of triangles argument");
    nt = di[0];
    nt2 = nt * 2;
    dnei = mxGetPr(prhs[1]);

    /* create internal neighbors list */
    nei = (int *) mxCalloc(MAXNRNE2 * nc, sizeof(int));
    if (nei == NULL)
        mexErrMsgTxt("Error allocating memory to store neighbors list.");
    nnei = (unsigned char *) mxCalloc(nc, sizeof(unsigned char));
    if (nnei == NULL)
        mexErrMsgTxt("Error allocating memory to store number of neighbors list.");
    nneic = (unsigned char *) mxCalloc(nc, sizeof(unsigned char));
    if (nneic == NULL)
        mexErrMsgTxt("Error allocating memory to store number of neighbors list.");

    /* create output argument */
    cdi[0] = nc;
    plhs[0] = mxCreateCellArray(2, cdi);
    if (plhs[0] == NULL)
        mexErrMsgTxt("Error allocating memory for cell array.");

    /* iterate over triangles */
    for (cc = 0; cc < nt; ++cc) {

        /* get three vertix numbers */
        dn1 = dnei[cc];
        dn2 = dnei[cc + nt];
        dn3 = dnei[cc + nt2];

        /* bail out if invalid */
        if (mxIsInf(dn1) || mxIsInf(dn2) || mxIsInf(dn3) ||
            mxIsNaN(dn1) || mxIsNaN(dn2) || mxIsNaN(dn3) ||
            dn1 < 1 || dn1 > nc || dn2 < 1 || dn2 > nc || dn3 < 1 || dn3 > nc)
            mexErrMsgTxt("Bad triangles argument.");

        /* get integer version */
        n1 = (int) (dn1 - 1);
        n2 = (int) (dn2 - 1);
        n3 = (int) (dn3 - 1);

        /* get offset to store neighbors of first vertex */
        sc = MAXNRNE2 * n1 + 2 * (nnei[n1]);
        ++nnei[n1];

        /* store neighbors */
        nei[sc++] = n2 + 1;
        nei[sc] = n3 + 1;

        /* repeat for second vertex */
        sc = MAXNRNE2 * n2 + 2 * (nnei[n2]);
        ++nnei[n2];

        nei[sc++] = n3 + 1;
        nei[sc] = n1 + 1;

        /* and third vertex also */
        sc = MAXNRNE2 * n3 + 2 * (nnei[n3]);
        ++nnei[n3];

        nei[sc++] = n1 + 1;
        nei[sc] = n2 + 1;
    }

    /* fill output array */
    cnei = plhs[0];
    cdi[0] = 1;
    ndi[1] = 1;
    for (cc = 0; cc < nc; ++cc) {

        /* set value */
        tnnei = nnei[cc];
        if (tnnei == 1)
            tnnei = 2;

        /* get offset into neighbors list array */
        neip = &nei[MAXNRNE2 * cc];

        /* get first neighbor */
        n1 = *neip;

        /* store into dnei buffer */
        *dneib = (double) n1;

        /* initialize counter to 1 */
        cinei = 1;

        /* get second neighbor and repeat until same as first */
        lcnt=0;
        for (n2 = neip[1]; n2 != n1; ) {
            if (lcnt++>1000)
              break;

            /* store next neighbor */
            dneib[cinei] = (double) n2;

            /* find next pair */
            for (cnnei = 0; cnnei < MAXNRNE2; cnnei += 2) {
                if (neip[cnnei] == n2) {
                    n2 = neip[++cnnei];
                    ++cinei;
                    break;
                }
            }
        }

        if (cinei != tnnei) {
            if (trackbn)
                bnl[bnc++] = cc + 1;
            sprintf(vstr, "Invalid neighborhood for vertex %d; cutting back.", cc + 1);
            mexWarnMsgTxt(vstr);
            tnnei = cinei;
        }

        /* create number of neighbors 1x1 array */
        onnei = mxCreateNumericArray(2, ndi, mxDOUBLE_CLASS, mxREAL);
        if (onnei == NULL)
            mexErrMsgTxt("Cannot allocate number of neighbors scalar.");
        *(mxGetPr(onnei)) = (double) tnnei;
        nnei[cc] = tnnei;

        /* create matching output array and copy */
        cdi[1] = tnnei;
        onei = mxCreateNumericArray(2, cdi, mxDOUBLE_CLASS, mxREAL);
        if (onei == NULL)
            mexErrMsgTxt("Cannot allocate list of neighbors array.");
        dnei = (double *) mxGetPr(onei);
        if (dnei == NULL)
            mexErrMsgTxt("Error getting double pointer to store neighbors.");
        for (cnnei = 0; cnnei < tnnei; ++cnnei)
            *dnei++ = dneib[cnnei];

        /* set to cell */
        mxSetCell(cnei, cc, onnei);
        mxSetCell(cnei, cc + nc, onei);
    }

}
