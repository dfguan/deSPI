/*
 * =====================================================================================
 *
 *       Filename:  kbtree.c
 *
 *    Description:  test kbtree
 *
 *        Version:  1.0
 *        Created:  04/26/2017 06:00:27 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#include <stdio.h>
#include "kbtree.h"
#include <stdint.h>
typedef struct {
	uint64_t key;
	int value;
	int hahah;
} intmap_t;

typedef struct {
	char *key;
	int value;
} strmap_t;

#define __intcmp(a, b) (((a).key > (b).key) - ((a).key < (b).key))
#define __strcmp(a, b) strcmp((a).key, (b).key)

#include "kbtree.h"
KBTREE_INIT(uint64_t, intmap_t, __intcmp)
/*KBTREE_INIT(str, strmap_t, __strcmp)*/

int test_int(int N, const unsigned *data)
{
	int i, ret;
	intmap_t *p, d;
	kbtree_t(uint64_t) *h;
	kbitr_t itr;
	h = kb_init(uint64_t, KB_DEFAULT_SIZE);
	for (i = 0; i < N; ++i) {
		d.key = data[i]; d.value = i;
		if (kb_getp(uint64_t, h, &d) == 0) kb_putp(uint64_t, h, &d);
		else kb_delp(uint64_t, h, &d);
		d.key = data[i], d.value = 1;
		p = kb_getp(uint64_t, h, &d); // kb_get() also works
		// IMPORTANT: put() only works if key is absent
        	if (!p) kb_putp(uint64_t, h, &d);
        	else ++p->value;
	
	}
    // ordered tree traversal
	kb_itr_first(uint64_t, h, &itr); // get an iterator pointing to the first
	for (; kb_itr_valid(&itr); kb_itr_next(uint64_t, h, &itr)) { // move on
		p = &kb_itr_key(intmap_t, &itr);
		printf("%d\t%lu\n", p->value, p->key);
	}		


	kb_destroy(int, h);
	return 0;
}
 
