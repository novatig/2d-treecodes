#include <cassert>
#include <cmath>

#include <algorithm>
#include <limits>
#include <vector>
#include <map>
#include <complex>
#include "treecode.h"

#define LEAF_MAXCOUNT 100
#define LMAX 15

#ifndef ORDER
#define ORDER 2
#endif

using namespace std;

namespace TreeCodeDiego
{
    const realtype eps = 10 * numeric_limits<realtype>::epsilon();

    realtype ext, xmin, ymin, thetasquared;

    vector<int> keys;
    vector<realtype> data[3];

    struct Node
    {
	int x, y, l, s, e;
	bool leaf;
	realtype xcom, ycom, r;
	
	complex<realtype> expansions[ORDER + 1];
    };

    map<int, Node> tree;

    void upward(Node& node)
    {    
	const int s = node.s, e = node.e;
    
	for(int i = s; i < e; ++i)
	{
	    const complex<realtype> rp = complex<realtype>(data[0][i] - node.xcom, data[1][i] - node.ycom);
	
	    complex<realtype> prod = rp;
	
	    for (int n = 1; n <= ORDER; ++n)
	    {
		node.expansions[n] -= prod * (data[2][i] / n); 

		prod *= rp;
	    }
	
	    node.expansions[0] += data[2][i];
	} 
    }

    constexpr unsigned int factorial(const int n)
    {
	return n <= 1 ? 1 : (n * factorial(n - 1));
    }

    constexpr unsigned int binomial(const int n, const int k)
    {
	return factorial(n) / (factorial(n - k) * factorial(k));
    }

    void upward(const Node src, Node& dst)
    {
	complex<realtype> rb(src.xcom - dst.xcom, src.ycom - dst.ycom);
    
	for (int j = 1; j <= ORDER; ++j)
	{
	    complex<realtype> csum, prod(1, 0);
	
	    for (int k = j; k >= 1; --k)
	    {
		csum += prod * src.expansions[k] * (realtype)binomial(j - 1, k - 1);

		prod *= rb;
	    }
	
	    csum -= prod * src.expansions[0] / (realtype)j;
	
	    dst.expansions[j] += csum;
	}

	dst.expansions[0] += src.expansions[0]; 
    }

    int nodeid(const int x, const int y, const int l)
    {
	return ((1 << 2 * l) - 1) / 3 + (x + (1 << l) * y);
    }

    Node build(const int x, const int y, const int l, const int s, const int e, const int mask)
    {
	const double h = ext / (1 << l);
	const double x0 = xmin + h * x, y0 = ymin + h * y;
	
	assert(x < (1 << l) && y < (1 << l) && x >= 0 && y >= 0);
    
#ifndef NDEBUG	
	for(int i = s; i < e; ++i)
	    assert(data[0][i] >= x0 && data[0][i] < x0 + h && data[1][i] >= y0 && data[1][i] < y0 + h);
#endif

	Node node = {x, y, l, s, e, e - s <= LEAF_MAXCOUNT || l + 1 > LMAX};
    
	{	
	    realtype xsum = 0, ysum = 0, weight = 0, r = 0;
	
	    for(int i = s; i < e; ++i)
		weight += fabs(data[2][i]);

	    for(int i = s; i < e; ++i)
		xsum += data[0][i] * fabs(data[2][i]);

	    for(int i = s; i < e; ++i)
		ysum += data[1][i] * fabs(data[2][i]);

	    if (weight != 0 && e - s > 0)
	    {
		node.xcom = xsum / weight;
		node.ycom = ysum / weight;
	    }
	    else
	    {
		node.xcom = x0 + h / 2;
		node.ycom = y0 + h / 2;
	    }
	    
	    for(int i = s; i < e; ++i)
		r = max(r, pow(data[0][i] - node.xcom, (realtype)2) + pow(data[1][i] - node.ycom, (realtype)2));
	
	    r = sqrt(r);
	    node.r = r;

	    assert(node.xcom >= x0 && node.xcom < x0 + h && node.ycom >= y0 && node.ycom < y0 + h);
	    assert(r < 1.5 * h);
	}
    
	if (node.leaf)
	    upward(node);
	else
	{
	    const vector<int>::const_iterator itbegins = keys.begin();
	    
	    for(int c = 0; c < 4; ++c)
	    {
		const int shift = 2 * (LMAX - l - 1);
	    
		const int key1 = mask | (c << shift);
		const int key2 = key1 + (1 << shift) - 1;

		const int indexmin = lower_bound(itbegins + s, itbegins + e, key1) - itbegins;
		const int indexsup = upper_bound(itbegins + s, itbegins + e, key2) - itbegins;

		Node child = build((x << 1) + (c & 1), (y << 1) + (c >> 1), l + 1, indexmin, indexsup, key1);

		upward(child, node);
	    }
	}
   
	const int entry = nodeid(x, y, l);
	assert(tree.find(entry) == tree.end());
    
	return tree[entry] = node;
    }

    realtype evaluate(const realtype xt, const realtype yt, const Node node)
    {
	const realtype r2 = pow(xt - node.xcom, 2) + pow(yt - node.ycom, 2);

	if (4 * node.r * node.r < thetasquared * r2)
	{
	    complex<realtype> z(xt - node.xcom, yt - node.ycom);
	    complex<realtype> s = node.expansions[0] * log(z);
	    complex<realtype> prod = complex<realtype>(1,0) / z;

	    for(int n = 1; n <= ORDER; ++n)
	    {
		s += prod * node.expansions[n];
		    
		prod /= z;
	    }

	    return s.real();
	}
	else
	{
	    realtype s = 0;

	    if (node.leaf)
		for(int i = node.s; i < node.e; ++i)
		{
		    const realtype xr = xt - data[0][i];
		    const realtype yr = yt - data[1][i];
		    s += log(sqrt(xr * xr + yr * yr + eps)) * data[2][i];
		}
	    else
		for(int c = 0; c < 4; ++c)
		{
		    const int entry = nodeid((node.x << 1) + (c & 1), (node.y << 1) + (c >> 1), node.l + 1);
		    assert(tree.find(entry) != tree.end());

		    s += evaluate(xt, yt, tree[entry]);
		}

	    return s;
	}
    }
}

using namespace TreeCodeDiego;

void treecode_potential(const realtype theta,
			const realtype * const xsrc, const realtype * const ysrc, const realtype * const vsrc, const int nsrc, 
			const realtype * const xdst, const realtype * const ydst, const int ndst, realtype * const vdst)
{    
    keys.resize(nsrc);
    
    for(int c = 0; c < 3; ++c)
	data[c].resize(nsrc);

    thetasquared = theta * theta;
    
    xmin = *min_element(xsrc, xsrc + nsrc);
    ymin = *min_element(ysrc, ysrc + nsrc);
    
    const realtype ext0 = (*max_element(xsrc, xsrc + nsrc) - xmin);
    const realtype ext1 = (*max_element(ysrc, ysrc + nsrc) - ymin);
    
    ext = max(ext0, ext1) * (1 + 2 * eps);
    xmin -= eps * ext;
    ymin -= eps * ext;

    vector< pair<int, int> > kv(nsrc);
    
    for(int i = 0; i < nsrc; ++i)
    {
	int x = floor((xsrc[i] - xmin) / ext * (1 << LMAX));
	int y = floor((ysrc[i] - ymin) / ext * (1 << LMAX));
	
	assert(x >= 0 && y >= 0);
	assert(x < (1 << LMAX) && y < (1 << LMAX));

	x = (x | (x << 8)) & 0x00FF00FF;
	x = (x | (x << 4)) & 0x0F0F0F0F;
	x = (x | (x << 2)) & 0x33333333;
	x = (x | (x << 1)) & 0x55555555;
	
	y = (y | (y << 8)) & 0x00FF00FF;
	y = (y | (y << 4)) & 0x0F0F0F0F;
	y = (y | (y << 2)) & 0x33333333;
	y = (y | (y << 1)) & 0x55555555;

	const int key = x | (y << 1);

	kv[i].first = key;
	kv[i].second = i;
    }

    sort(kv.begin(), kv.end());
    
    for(int i = 0; i < nsrc; ++i)
    {
	keys[i] = kv[i].first;

	const int entry = kv[i].second;
	assert(entry >= 0 && entry < nsrc);
	
	data[0][i] = xsrc[entry];
	data[1][i] = ysrc[entry];
	data[2][i] = vsrc[entry];
    }

    kv.clear();
    tree.clear();

    Node root = build(0, 0, 0, 0, nsrc, 0);

    for(int i = 0; i < ndst; ++i)
	vdst[i] = evaluate(xdst[i], ydst[i], root);
    
    for(int c = 0; c < 3; ++c)
	data[c].clear();
}

