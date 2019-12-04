//=================================================================================
//  adaptels.cpp
//=================================================================================
/*
* Copyright (c) 2019, Radhakrishna Achanta, Ecole Polytechnique Federale de Lausanne (EPFL)
* All rights reserved.
* Redistribution and use in source and binary forms, with or without
* modification, are permitted ONLY FOR NON-COMMERCIAL PURPOSES provided that
* the following conditions are met:
*
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*        notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of EPFL nor the author may be used to endorse or
*       promote products derived from this software without specific prior
*       written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE AUTHOR "AS IS" AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL EPFL OR THE AUTHOR BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*
* FOR ANY COMMERCIAL USE, PLEASE CONTACT THE AUTHOR AT: radhakrishna.achanta@epfl.ch
*/

#include <stdlib.h>
#include <math.h>
// #include <stdbool.h>
#include "adaptels.h"

void push (HEAP *h, const int ind, const int x, const int y, const double dist)
{
    if (h->len + 1 >= h->size)
    {
        h->size = h->size ? h->size * 2 : 4; // if h->size is greater than 0 mulptiply by 2 else set it to 4 (why?)
        h->nodes = (NODE *)realloc(h->nodes, h->size * sizeof (NODE));
    }
    int i = h->len + 1;
    int j = i / 2;
    while (i > 1 && h->nodes[j].d > dist)
    {
        h->nodes[i] = h->nodes[j];
        i = j;
        j = j / 2;
    }
    h->nodes[i].i = ind;
    h->nodes[i].x = x;
    h->nodes[i].y = y;
    h->nodes[i].d = dist;
    h->len++;
}
 
void pop (HEAP *h, NODE* pnode)
{
    int i, j, k;
    if (!h->len)
    {
        // return NULL;
        return;
    }
    //int i = h->nodes[1].i;
    pnode->i = h->nodes[1].i;
    pnode->x = h->nodes[1].x;
    pnode->y = h->nodes[1].y;
    pnode->d = h->nodes[1].d;
 
    h->nodes[1] = h->nodes[h->len];
 
    h->len--;
 
    i = 1;
    while (i!=h->len+1)
    {
        k = h->len+1;
        j = 2 * i;
        if (j <= h->len && h->nodes[j].d < h->nodes[k].d)
        {
            k = j;
        }
        if (j + 1 <= h->len && h->nodes[j + 1].d < h->nodes[k].d)
        {
            k = j + 1;
        }
        h->nodes[i] = h->nodes[k];
        i = k;
    }
    // return i;
}

// void rgbtolab(int* rin, int* gin, int* bin, int sz, double* lvec, double* avec, double* bvec)
void rgbtolab(double* rin, double* gin, double* bin, int sz, double* lvec, double* avec, double* bvec)
{
    int i;
    double sR, sG, sB;
    double R,G,B;
    double X,Y,Z;
    double r, g, b;
    const double epsilon = 0.008856;	//actual CIE standard
    const double kappa   = 903.3;		//actual CIE standard
    
    const double Xr = 0.950456;	//reference white
    const double Yr = 1.0;		//reference white
    const double Zr = 1.088754;	//reference white
    double xr,yr,zr;
    double fx, fy, fz;
    double lval,aval,bval;
    
    for(i = 0; i < sz; i++)
    {
        sR = rin[i]; sG = gin[i]; sB = bin[i];
        R = sR/255.0;
        G = sG/255.0;
        B = sB/255.0;
        
        if(R <= 0.04045)	r = R/12.92;
        else				r = pow((R+0.055)/1.055,2.4);
        if(G <= 0.04045)	g = G/12.92;
        else				g = pow((G+0.055)/1.055,2.4);
        if(B <= 0.04045)	b = B/12.92;
        else				b = pow((B+0.055)/1.055,2.4);
        
        X = r*0.4124564 + g*0.3575761 + b*0.1804375;
        Y = r*0.2126729 + g*0.7151522 + b*0.0721750;
        Z = r*0.0193339 + g*0.1191920 + b*0.9503041;
        
        //------------------------
        // XYZ to LAB conversion
        //------------------------
        xr = X/Xr;
        yr = Y/Yr;
        zr = Z/Zr;
        
        if(xr > epsilon)	fx = pow(xr, 1.0/3.0);
        else				fx = (kappa*xr + 16.0)/116.0;
        if(yr > epsilon)	fy = pow(yr, 1.0/3.0);
        else				fy = (kappa*yr + 16.0)/116.0;
        if(zr > epsilon)	fz = pow(zr, 1.0/3.0);
        else				fz = (kappa*zr + 16.0)/116.0;
        
        lval = 116.0*fy-16.0;
        aval = 500.0*(fx-fy);
        bval = 200.0*(fy-fz);
        
        lvec[i] = lval; avec[i] = aval; bvec[i] = bval;
    }
}

//===========================================================================
/// createAdaptels
//===========================================================================
void createAdaptels(
             // double*		lv,
             // double*		av,
             // double*		bv,
            double**                    chans,
            const int                   nchans,
             const int					width,
             const int					height,
             int*				labels,
             int*						numklabels,
             const double               threshold)
{
    // Settings
    const int CONNECTIVITY = 8;
    // const bool useAccumulatedAverage = true;//setting this to false gives poor curves
    // const double SMALLCONST = 0;//either 1 or 0
    const double THRESH = threshold;// for accumulated distance along the geodesic paths
    
    // Initializations
    const int w = width;
    const int h = height;
    const int sz = w*h;
    const int dx8[8] = {-1,  0, 1, 0, -1,  1, 1, -1};//for 4 or 8 connectivity
    const int dy8[8] = { 0, -1, 0, 1, -1, -1, 1,  1};//for 4 or 8 connectivity
    const int dn8[8] = {-1, -w, 1, w, -1-w,1-w,1+w,-1+w};
    // const double xy8[8] = {1,1,1,1,2,2,2,2};
    
    // struct NODE
    // {
    //     int x,y,i;
    //     double d;
    //     bool operator()(const NODE& one, const NODE& two)
    //     {
    //         return one.d > two.d;
    //     }
    // };
    
    // vector<double> distvec(sz,0);//this is the distance map that store the accumulated distances
    double* distvec = (double*)malloc(sizeof(double)*sz);
    for(int i = 0; i < sz; i++) labels[i] = -1;
    int numLabels = 0;
    
    // vector<int> xvec(sz),yvec(sz),indexvec(sz);//the candidates for starting a new segment
    int* xvec = (int*)malloc(sizeof(int)*sz);
    int* yvec = (int*)malloc(sizeof(int)*sz);
    int* indexvec = (int*)malloc(sizeof(int)*sz);
    //vector<bool> istaken(sz,false);
    xvec[0] = width>>1; yvec[0] = height>>1;
    indexvec[0] = yvec[0]*width+xvec[0]; int indexcount = 1;
    
    // priority_queue<NODE, vector<NODE>, NODE> pq;
    // NODE node,temp;

    HEAP* pheap = (HEAP *)calloc(1, sizeof (HEAP));//create a minheap
    NODE* pnode = (NODE*)calloc(1, sizeof(NODE));
    int i, n, p, c, oindex, segsz, heapsize;
    int xx, yy, ii;
    // double lsum,asum,bsum;
    // double ldiff,adiff,bdiff;
    double* csum = (double*)malloc(sizeof(double)*nchans);
    double cdiff = 0;
    double colordist,distsum;
    for(p = 0; p < indexcount; p++)
    {
        oindex = indexvec[p];
        if(labels[oindex] < 0)
        {
            // node.d = 0; node.x = xvec[p]; node.y = yvec[p]; node.i = oindex;
            // pq.push(node);
            push(pheap,oindex,xvec[p],yvec[p],0);
            distvec[oindex] = 0;
            labels[oindex] = numLabels;
            
            // lsum = lv[oindex]; asum = av[oindex]; bsum = bv[oindex];
            for(c = 0; c < nchans; c++)
            {
                csum[c] = chans[c][oindex];
            }
            segsz = 1;
            heapsize = 1;
            while(heapsize)//start a new adaptel
            {
                // node = pq.top(); pq.pop();heapsize--;
                pop(pheap,pnode);heapsize--;
                i = pnode->i;
                
                for(n = 0; n < CONNECTIVITY; n++)
                {
                    xx = pnode->x + dx8[n];
                    yy = pnode->y + dy8[n];
                    if(!(xx < 0 || xx >= w || yy < 0 || yy >= h))
                    {
                        ii = pnode->i + dn8[n];
                        if(distvec[i] < THRESH)
                        {
                            if(labels[ii] != labels[i])//to prevent relabeling pixels with the same labels
                            {
                                // ldiff = lsum - lv[ii]*segsz; adiff = asum - av[ii]*segsz; bdiff = bsum - bv[ii]*segsz;
                                // dist = sqrt(ldiff*ldiff + adiff*adiff + bdiff*bdiff)/segsz;
                                colordist = 0;
                                for(c = 0; c < nchans; c++)
                                {
                                    cdiff = csum[c]-(chans[c][ii]*segsz);
                                    colordist += (cdiff*cdiff);
                                }
                                colordist = sqrt(colordist)/segsz;
                                distsum = distvec[i] + colordist;

                                if(distsum < distvec[ii] || labels[ii] < 0)
                                {
                                    distvec[ii] = distsum;
                                    labels[ii] = labels[i];
                                    // lsum += lv[ii];asum += av[ii];bsum += bv[ii];segsz++;//accumulate color values
                                    for(c = 0; c < nchans; c++)
                                    {
                                        csum[c] += chans[c][ii];
                                    }
                                    segsz++;
                                    // temp.x = xx; temp.y = yy; temp.i = ii;temp.d = distvec[ii];
                                    // pq.push(temp);
                                    // heapsize++;
                                    push(pheap,ii,xx,yy,distvec[ii]); heapsize++;
                                }
                            }
                        }
                        else//put pixel index into the seed buffer
                        {
                            if(labels[ii] < 0)//i.e. take only unlabeled pixels as new starting points
                            {
                                indexvec[indexcount] = ii;
                                xvec[indexcount] = xx;
                                yvec[indexcount] = yy;
                                indexcount++;
                            }
                        }
                    }
                }//for n neighbours
            }//while loop ends
            numLabels++;
        }
    }
    *numklabels = numLabels;

    free(distvec);
    free(xvec);
    free(yvec);
    free(indexvec);
    free(csum);
    free(pnode);
    free(pheap->nodes);
    free(pheap);
}

//===========================================================================
/// Adaptels_main
///
/// The main function
//===========================================================================
void Adaptels_main(double* img, const int width, const int height,
                const int nchannels, const double T,
                const int doRGBtoLAB, int* klabels, int* numlabels)
{
    int sz = width*height;
    double** channels = (double**)malloc(sizeof(double*)*nchannels);
    for(int c = 0; c < nchannels; c++)
    {
        channels[c] = img + c*sz;
    }
    //---------------------------
    // Perform color conversion
    //---------------------------
    if(doRGBtoLAB && nchannels==3)
    {
        rgbtolab(channels[0],channels[1],channels[2],sz,channels[0],channels[1],channels[2]);
    }
    //---------------------------
    // Create adaptels
    //---------------------------
    createAdaptels(channels,nchannels,width,height,klabels,numlabels,T);

    free(channels);
}




