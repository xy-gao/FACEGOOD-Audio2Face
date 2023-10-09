// LPC.cpp : 定义 DLL 应用程序的导出函数。
//
/*!
 * \file
 * \brief Implementations of linear prediction functions, and conversion
 * between common representations of linear predictive parameters
 * \author Thomas Eriksson
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#include "lpcfunc.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LPC.h"

extern "C" {
	Vec *createVec(int size)
	{
		if (size <= 0)
			return NULL;
		Vec *vec = (Vec *)malloc(sizeof(Vec));
		if (vec)
		{
			vec->p = (double *)malloc(sizeof(double) * size);
			if (!vec->p)
			{
				free(vec);
				vec = NULL;
			}
			else
			{
				vec->size = size;
				memset(vec->p, 0, vec->size * sizeof(double));
			}
		}
		return vec;
	}
	void clearVec(Vec *vec)
	{
		memset(vec->p, 0, sizeof(double) * vec->size);
	}


	void deleteVec(Vec *vec)
	{
		free(vec->p);
		vec->p = NULL;
		free(vec);
	}

	Vec* autocorr(const Vec *x, int order)
	{
		if (order < 0) order = x->size;

		Vec *R = createVec(order + 1);
		if (!R)
			return NULL;
		double sum;
		int i, j;

		for (i = 0; i < order + 1; i++) {
			sum = 0;
			for (j = 0; j < x->size - i; j++) {
				sum += x->p[j] * x->p[j + i];
			}
			R->p[i] = sum;
		}
		return R;
	}

	Vec* levinson(const Vec *R2, int order)
	{
		Vec *R = (Vec*)R2;
		R->p[0] = R->p[0] * (1. + 1.e-9);

		if (order < 0) order = R->size - 1;
		double k, alfa, s;
		double *any = (double *)malloc(sizeof(double) * (order + 1));
		double *a = (double *)malloc(sizeof(double) * (order + 1));
		int j, m;
		Vec *out = createVec(order + 1);

		a[0] = 1;
		alfa = R->p[0];
		if (alfa <= 0) {
			out->p[0] = 1;
			return out;
		}

		for (m = 1; m <= order; m++) {
			s = 0;
			for (j = 1; j < m; j++) {
				s = s + a[j] * R->p[m - j];
			}

			k = -(R->p[m] + s) / alfa;
			if (fabs(k) >= 1.0) {
				printf("levinson : panic! abs(k)>=1, order &d. Aborting...\n");
				for (j = m; j <= order; j++) {
					a[j] = 0;
				}
				break;
			}
			for (j = 1; j < m; j++) {
				any[j] = a[j] + k * a[m - j];
			}
			for (j = 1; j < m; j++) {
				a[j] = any[j];
			}
			a[m] = k;
			alfa = alfa * (1 - k * k);
		}
		for (j = 0; j < out->size; j++) {
			out->p[j] = a[j];
		}
		free(any);
		free(a);
		return out;
	}

	//void LPC(const double *x, int size, int order, double *out)
	//{
	//	Vec *vec = createVec(size);
	//	memcpy(vec, x, size * sizeof(double));

	//	Vec *temp = levinson(autocorr(vec, order - 1), order - 1);
	//	memcpy(out, temp->p, temp->size);
	//	free(vec->p);
	//	free(vec);
	//	free(temp->p);
	//	free(temp);
	//	//return d;
	//}
 int LPC(double *in, int size, int order, double *out)
	{
		//printf("LPC(float *x, int size, int order, float *out)");
		Vec *vec = createVec(size);
		for (int i = 0; i < vec->size; i++)
			vec->p[i] = in[i];
		Vec *temp = levinson(autocorr(vec, order - 1), order - 1);
		for (int i = 0; i < temp->size; i++)
			out[i] = temp->p[i];

		free(vec->p);
		free(vec);
		free(temp->p);
		free(temp);
		return order;
		//return d;
	}

	Vec *poly2ac(const Vec *poly)
	{
		Vec *a = (Vec*)poly;
		int order = a->size - 1;
		double alfa, s, *any = (double *)malloc((order + 1) * sizeof(double));
		int j, m;
		Vec *r = createVec(order + 1);
		Vec *k = poly2rc(a);

		if (a->p[0] != 1)
			printf("poly2ac : not an lpc filter\n");
		r->p[0] = 1;
		alfa = 1;
		for (m = 1; m <= order; m++) {
			s = 0;
			for (j = 1; j < m; j++) {
				s = s + a->p[j] * r->p[m - j];
			}
			r->p[m] = -s - alfa * k->p[m - 1];
			for (j = 1; j < m; j++) {
				any[j] = a->p[j] + k->p[m - 1] * a->p[m - j];
			}
			for (j = 1; j < m; j++) {
				a->p[j] = any[j];
			}
			a[m] = k[m - 1];
			alfa = alfa * (1 - sqr(k->p[m - 1]));
		}
		free(any);
		return r;
	}

	Vec *poly2rc(const Vec *a)
	{
		// a is [1 xx xx xx], a.size()=order+1
		int   m, i;
		int    order = a->size - 1;
		Vec *k = createVec(order);
		Vec *any = createVec(order + 1), *aold = createVec(a->size);

		for (m = order - 1; m > 0; m--) {
			k->p[m] = aold->p[m + 1];
			if (fabs(k->p[m]) > 1) k->p[m] = 1.0 / k->p[m];
			for (i = 0; i < m; i++) {
				any->p[i + 1] = (aold->p[i + 1] - aold->p[m - i] * k->p[m]) / (1 - k->p[m] * k->p[m]);
			}
			aold = any;
		}
		k->p[0] = any->p[1];
		if (fabs(k->p[0]) > 1) k->p[0] = 1.0 / k->p[0];
		return k;
	}

	Vec *rc2poly(const Vec *k)
	{
		int  m, i;
		Vec *a = createVec(k->size + 1), *any = createVec(k->size + 1);

		a->p[0] = 1;
		any->p[0] = 1;
		a->p[1] = k->p[0];
		for (m = 1; m < k->size; m++) {
			any[m + 1] = k[m];
			for (i = 0; i < m; i++) {
				any->p[i + 1] = a->p[i + 1] + a->p[m - i] * k->p[m];
			}
			a = any;
		}
		return a;
	}

	Vec *rc2lar(const Vec *k)
	{
		short m;
		Vec *LAR = createVec(k->size);

		for (m = 0; m < k->size; m++) {
			LAR->p[m] = log((1 + k->p[m]) / (1 - k->p[m]));
		}
		return LAR;
	}

	Vec *lar2rc(const Vec *LAR)
	{
		short m;
		Vec *k = createVec(LAR->size);

		for (m = 0; m < LAR->size; m++) {
			k->p[m] = (exp(LAR->p[m]) - 1) / (exp(LAR->p[m]) + 1);
		}
		return k;
	}

	double FNevChebP_double(double  x, const double c[], int n)
	{
		int i;
		double b0 = 0.0, b1 = 0.0, b2 = 0.0;

		for (i = n - 1; i >= 0; --i) {
			b2 = b1;
			b1 = b0;
			b0 = 2.0 * x * b1 - b2 + c[i];
		}
		return (0.5 * (b0 - b2 + c[0]));
	}

	double FNevChebP(double  x, const double c[], int n)
	{
		int i;
		double b0 = 0.0, b1 = 0.0, b2 = 0.0;

		for (i = n - 1; i >= 0; --i) {
			b2 = b1;
			b1 = b0;
			b0 = 2.0 * x * b1 - b2 + c[i];
		}
		return (0.5 * (b0 - b2 + c[0]));
	}

	Vec *poly2lsf(const Vec *pc)
	{
		int np = pc->size - 1;
		Vec *lsf = createVec(np);

		Vec *fa = createVec((np + 1) / 2 + 1), *fb = createVec((np + 1) / 2 + 1);
		Vec *ta = createVec((np + 1) / 2 + 1), *tb = createVec((np + 1) / 2 + 1);
		double *t;
		double xlow, xmid, xhigh;
		double ylow, ymid, yhigh;
		double xroot;
		double dx;
		int i, j, nf;
		int odd;
		int na, nb, n;
		double ss, aa;
		double DW = (0.02 * dPi);
		int  NBIS = 4;

		odd = (np % 2 != 0);
		if (odd) {
			nb = (np + 1) / 2;
			na = nb + 1;
		}
		else {
			nb = np / 2 + 1;
			na = nb;
		}

		fa->p[0] = 1.0;
		for (i = 1, j = np; i < na; ++i, --j)
			fa->p[i] = pc->p[i] + pc->p[j];

		fb->p[0] = 1.0;
		for (i = 1, j = np; i < nb; ++i, --j)
			fb->p[i] = pc->p[i] - pc->p[j];

		if (odd) {
			for (i = 2; i < nb; ++i)
				fb->p[i] = fb->p[i] + fb->p[i - 2];
		}
		else {
			for (i = 1; i < na; ++i) {
				fa->p[i] = fa->p[i] - fa->p[i - 1];
				fb->p[i] = fb->p[i] + fb->p[i - 1];
			}
		}

		ta->p[0] = fa->p[na - 1];
		for (i = 1, j = na - 2; i < na; ++i, --j)
			ta->p[i] = 2.0 * fa->p[j];

		tb[0] = fb[nb - 1];
		for (i = 1, j = nb - 2; i < nb; ++i, --j)
			tb->p[i] = 2.0 * fb->p[j];

		nf = 0;
		t = ta->p;
		n = na;
		xroot = 2.0;
		xlow = 1.0;
		ylow = FNevChebP_double(xlow, t, n);


		ss = sin(DW);
		aa = 4.0 - 4.0 * cos(DW) - ss;
		while (xlow > -1.0 && nf < np) {
			xhigh = xlow;
			yhigh = ylow;
			dx = aa * xhigh * xhigh + ss;
			xlow = xhigh - dx;
			if (xlow < -1.0)
				xlow = -1.0;
			ylow = FNevChebP_double(xlow, t, n);
			if (ylow * yhigh <= 0.0) {
				dx = xhigh - xlow;
				for (i = 1; i <= NBIS; ++i) {
					dx = 0.5 * dx;
					xmid = xlow + dx;
					ymid = FNevChebP_double(xmid, t, n);
					if (ylow * ymid <= 0.0) {
						yhigh = ymid;
						xhigh = xmid;
					}
					else {
						ylow = ymid;
						xlow = xmid;
					}
				}
				if (yhigh != ylow)
					xmid = xlow + dx * ylow / (ylow - yhigh);
				else
					xmid = xlow + dx;
				lsf->p[nf] = acos((double)xmid);
				++nf;
				if (xmid >= xroot) {
					xmid = xlow - dx;
				}
				xroot = xmid;
				if (t == ta->p) {
					t = tb->p;
					n = nb;
				}
				else {
					t = ta->p;
					n = na;
				}
				xlow = xmid;
				ylow = FNevChebP_double(xlow, t, n);
			}
		}
		if (nf != np) {
			printf("poly2lsf: WARNING: failed to find all lsfs\n");
		}
		return lsf;
	}

	Vec *lsf2poly(const Vec *f)
	{
		int m = f->size;
		Vec  *pc = createVec(m + 1);
		double c1, c2, *a;
		Vec  *p = createVec(m + 1), *q = createVec(m + 1);
		int mq, n, i, nor;

		if (m % 2 != 0)
			printf("lsf2poly: THIS ROUTINE WORKS ONLY FOR EVEN m\n");
		pc->p[0] = 1.0;
		a = pc->p + 1;
		mq = m >> 1;
		for (i = 0; i <= m; i++) {
			q->p[i] = 0.;
			p->p[i] = 0.;
		}
		p->p[0] = q->p[0] = 1.;
		for (n = 1; n <= mq; n++) {
			nor = 2 * n;
			c1 = 2 * cos(f->p[nor - 1]);
			c2 = 2 * cos(f->p[nor - 2]);
			for (i = nor; i >= 2; i--) {
				q->p[i] += q->p[i - 2] - c1 * q->p[i - 1];
				p->p[i] += p->p[i - 2] - c2 * p->p[i - 1];
			}
			q->p[1] -= c1;
			p->p[1] -= c2;
		}
		a[0] = 0.5 * (p->p[1] + q->p[1]);
		for (i = 1, n = 2; i < m; i++, n++)
			a[i] = 0.5 * (p->p[i] + p->p[n] + q->p[n] - q->p[i]);

		return pc;
	}

	Vec *poly2cepstrum(const Vec *a)
	{
		Vec *c = createVec(a->size - 1);

		for (int n = 1; n <= c->size; n++) {
			c[n - 1] = a[n];
			for (int k = 1; k < n; k++) {
				c->p[n - 1] -= ((double)(k)) / n * a->p[n - k] * c->p[k - 1];
			}
		}
		return c;
	}

	Vec *poly2cepstrumVI(const Vec *a, int num)
	{
		if (num < a->size)
			printf("a2cepstrum : not allowed cepstrum length\n");
		Vec *c = createVec(num);
		int n;

		for (n = 1; n < a->size; n++) {
			c[n - 1] = a[n];
			for (int k = 1; k < n; k++) {
				c->p[n - 1] -= ((double)(k)) / n * a->p[n - k] * c->p[k - 1];
			}
		}
		for (n = a->size; n <= c->size; n++) {
			c->p[n - 1] = 0;
			for (int k = n - a->size + 1; k < n; k++) {
				c->p[n - 1] -= ((double)(k)) / n * a->p[n - k] * c->p[k - 1];
			}
		}
		return c;
	}

	Vec *cepstrum2poly(const Vec *c)
	{
		Vec *a = createVec(c->size + 1);

		a->p[0] = 1;
		for (int n = 1; n <= c->size; n++) {
			a[n] = c[n - 1];
			for (int k = 1; k < n; k++) {
				a->p[n] += ((double)(k)) / n * a->p[n - k] * c->p[k - 1];
			}
		}
		return a;
	}

	Vec *chirp(const Vec *a, double factor)
	{
		Vec *temp = createVec(a->size);
		int i;
		double   f = factor;

		if (a->p[0] != 1)
			printf("chirp : a[0] should be 1");
		temp->p[0] = a->p[0];
		for (i = 1; i < a->size; i++) {
			temp->p[i] = a->p[i] * f;
			f *= factor;
		}
		return temp;
	}

	Vec *schurrc(const Vec *R, int order)
	{
		if (order == -1) order = R->size - 1;

		Vec *k = createVec(order), *scratch = createVec(2 * order + 2);

		int m;
		int h;
		double ex;
		double *ep;
		double *en;

		ep = scratch->p;
		en = scratch->p + order + 1;

		m = 0;
		while (m < order) {
			m++;
			ep[m] = R->p[m];
			en[m] = R->p[m - 1];
		}
		if (en[1] < 1.0) en[1] = 1.0;
		h = -1;
		while (h < order) {
			h++;
			k->p[h] = -ep[h + 1] / en[1];
			en[1] = en[1] + k->p[h] * ep[h + 1];
			if (h == (order - 1)) {
				// cout << "k: " << k << endl ;
				return k;
			}
			ep[order] = ep[order] + k->p[h] * en[order - h];
			m = h + 1;
			while (m < (order - 1)) {
				m++;
				ex = ep[m] + k->p[h] * en[m - h];
				en[m - h] = en[m - h] + k->p[h] * ep[m];
				ep[m] = ex;
			}
		}
		return k;  // can never come here
	}

	Vec *lerouxguegenrc(const Vec *R, int order)
	{
		Vec *k = createVec(order);

		double  *r, *rny;
		int j, m;
		int M = order;

		r = (double *)malloc(sizeof(double) * (2 * M + 1));
		rny = (double *)malloc(sizeof(double) * (2 * M + 1));

		for (j = 0; j <= M; j++) {
			r[M - j] = r[M + j] = R->p[j];
		}
		for (m = 1; m <= M; m++) {
			k->p[m - 1] = -r[M + m] / r[M];
			for (j = -M; j <= M; j++) {
				rny[M + j] = r[M + j] + k->p[m - 1] * r[M + m - j];
			}
			for (j = -M; j <= M; j++) {
				r[M + j] = rny[M + j];
			}
		}
		free(r);
		free(rny);
		return k;
	}
}
//double sd(const vec &In1, const vec &In2)
//{
//  return std::sqrt(37.722339402*energy(poly2cepstrum(In1, 32) - poly2cepstrum(In2, 32)));
//}
//
//// highestfreq=1 gives entire band
//double sd(const vec &In1, const vec &In2, double highestfreq)
//{
//  vec Diff = sqr(abs(log10(filter_spectrum(In1, In2))));
//  double S = 0;
//  for (int i = 0;i < round(highestfreq*129);i++) {
//    S = S + Diff(i);
//  }
//  S = S * 100 / round(highestfreq * 129);
//  return std::sqrt(S);
//}

//! \endcond
