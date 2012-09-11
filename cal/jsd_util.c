#include "jsd_util.h"

int jsd_bsearch(const void *key, const void *base, size_t nmemb,
              size_t size, int(*compar)(const void *, const void *))
{
	int low = -1;
	int high = nmemb;
	int mid;
	int comp;

//printf("jsd_search\n");
	while (high-low > 1) 
	{
		mid = (high+low)/2;
		comp = compar(key, base+size*mid);

//printf("low=%i high=%i mid=%i comp=%i\n", low,high,mid,comp);
		if (comp > 0) low = mid;
		else if (comp < 0) high = mid;
		else return mid;
	}
//printf("not found\n");
	return -1;
}
