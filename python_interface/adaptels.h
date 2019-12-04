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

typedef struct
{
    int x,y,i;
    double d;
} NODE;

typedef struct
{
    NODE *nodes;
    int len;
    int size;
} HEAP;

void push (HEAP *h, const int ind, const int x, const int y, const double dist);
void pop (HEAP *h, NODE* pnode);
void rgbtolab(double* rin, double* gin, double* bin, int sz, double* lvec, double* avec, double* bvec);
void createAdaptels(
             // double*		lv,
             // double*		av,
             // double*		bv,
            double**                    chans,
            const int                   nchans,
             const int					width,
             const int					height,
             int*				labels,
             int*						numlabels,
             const double               threshold);
void Adaptels_main(double* img, const int width, const int height,
                const int nchannels, const double T,
                const int doRGBtoLAB, int* labels, int* numlabels);