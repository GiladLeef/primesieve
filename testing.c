#include <time.h>
#include <sys/time.h>
// Basic ~ 100 lines factorization tester, use test=1 or test=160 to execute the tests.

static inline void facMiniTests(facParams *m) {
	// Welcome to the testing feature, it uses 2 numbers + a computation sheet

	cint nums[2];
	for (int i = 0; i < 2; ++i)
		cint_init(&nums[i], 2048, 0);
	cint *n = &nums[0], *factor = n + 1;
	cint_sheet *sheet = cint_new_sheet(1 << 10);
	facParams params = {0}; params.silent = 1 ;
	//---------------

	unsigned sr = add_rand_seed(sheet);
	sr ^= time(0);

	qs_sm errorNumber = 0, nth, bits = m->testing >= 3 && m->testing <= 300 ? (qs_sm) m->testing : 1;
	qs_sm seconds = bits < 140 ? 30 : bits < 160 ? 60 : bits < 180 ? 120 : 180, trialMax;
	const char *secondsStr = bits < 140 ? "30 seconds" : bits < 160 ? "minute" : bits < 180 ? "2 minutes" : "3 minutes";
	if (m->testing >= 3)  printf("-- %3d-bit : your %s factorization test -- \n\n", bits, secondsStr);
	else printf("-- crescendo : your %s factorization test -- \n\n", secondsStr);

	double chronometer = 0, timeout = 1e6 * seconds, took;

	for ( nth = 1; chronometer < timeout && !errorNumber && nth <= 1000; ++nth) {
		if (m->testing < 3) bits = nth > 2 ? nth : 3 ;
		params.silent = bits < 150 ;
		params.qs_limit = bits ;

		if (bits < 40) trialMax = 0 ;
		else if (bits <= 64) trialMax = 1024;
		else {
			trialMax = 4669921;
			for (qs_sm i = 0; i < 250; trialMax >>= (bits < i) * 1, i += 30);
		}

		retry :
		cint_random_bits(n, bits), *n->mem |= trialMax != 0;
		if (cint_is_prime(sheet, n, 2))
			goto retry; // a prime number isn't submitted

		for (qs_sm nVal = 3; nVal < trialMax; nVal += 2)
			if (is_prime_4669921(nVal))
				if (simple_int_to_cint(factor, nVal), cint_remove(sheet, n, factor))
					goto retry; // a trial divisible number isn't submitted

		char *str = cint_to_string(n, 10);

		printf(params.silent ? "%2d. %s = " : "%2d. %s\n", nth, str);
		free(str);
		fflush(stdout);

		struct timeval tv;
		gettimeofday(&tv, 0), took = -(tv.tv_sec * 1e6 + tv.tv_usec);
		facCint **factors = cFactor(n, &params);
		assert(factors); // answer is not null
		gettimeofday(&tv, 0), took += tv.tv_sec * 1e6 + tv.tv_usec;
		chronometer += took ;

		sr += (unsigned) timeout;
		srand(sr);

		if (params.silent == 0) printf("    [ %.2fs ] ", took / 1e6);

		for (int i = 0; factors[i]; ++i) {
			str = cint_to_string(&factors[i]->cint, 10);
			// a power of the answered factor must be removable from n.
			const unsigned powers = cint_remove(sheet, n, &factors[i]->cint);
			switch (powers) {
				case 0 : errorNumber |= 2; printf(" [%s] ", str); break;
				case 1 : printf(factors[i + 1] ? "%s * " : "%s", str); break;
				default: printf(factors[i + 1] ? "(%s ^ %d) * " : "(%s ^ %d)", str, powers); break;
			}
			free(str);
			fflush(stdout);
			const int nBits = (int) cint_count_bits(&factors[i]->cint);
			if (nBits != factors[i]->bits) errorNumber |= 4;

			if (nBits > 1) {
				// the factor must be 'prime' with 3 Miller-Rabin iterations.
				const int prime = (int) cint_is_prime(sheet, &factors[i]->cint, 3);
				if (prime == 0) errorNumber |= 8;
			}

		}

		// number n must be equal to 1 after all factors removed.
		if (cint_count_bits(n) != 1) errorNumber |= 16;
		printf(params.silent ? "\n" : "\n\n");

		free(factors); // clears this answer.
	}

	if (errorNumber == 0 && m->testing >= 3) puts("Prime numbers wasn't submitted, the software has checked the answers.");
	if (errorNumber & 2) puts("number wasn't a multiple of the factor");
	if (errorNumber & 4) puts("bit count of the factor was wrong");
	if (errorNumber & 8) puts("factor of number wasn't 'prime'");
	if (errorNumber & 16) puts("number wasn't correctly factored");
	if (errorNumber) printf("qs rand seed was %u\n", params.qs_rand_seed);

	chronometer /= 1e6 ;
	printf("Thank you, technical chronometer displays    %6.3f s   .\n", chronometer);
	if (m->testing >= 3) printf("On average a %3d-bit factorization take      %6.3f s   .\n",  bits, chronometer / nth);

	// Clear the numbers + the computation sheet.
	cint_clear_sheet(sheet);
	for (int i = 0; i < 2; ++i)
		free(nums[i].mem);
}
