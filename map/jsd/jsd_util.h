#ifndef _JSD_UTIL_H
#define _JSD_UTIL_H

#include <stdlib.h>


/* binary search
 * Emulates the stdlib bsearch, but returns the index of the match.
 * The array must be sorted, and should be sorted by the given 
 * compar function.
 *
 * @param key the element to use with the compar function
 * @param base the base of the sorted array
 * @param nmemb the number of members in the array
 * @param size the size of the array elements
 * @param compar the comparison function
 * @return -1 if not found or index of a matching item
 */
int jsd_bsearch(const void *key, const void *base, size_t nmemb,
              size_t size, int(*compar)(const void *, const void *));

#endif //_JSD_UTIL_H

