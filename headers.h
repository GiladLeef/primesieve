#ifndef FAC_HEADERS
#define FAC_HEADERS

#include <errno.h>
#include <stddef.h>

// Quadratic sieve integers (it has sometimes been tested with signed integers too).
typedef uint32_t qs_sm; // small size, like a factor base prime number (32-bit)
typedef uint64_t qs_md; // medium size, like a factor base prime number squared (64-bit)
typedef int64_t qs_tmp; // signed type to perform intermediates computations.

// The factorization manager (calls the quadratic sieve)

typedef struct {
	cint cint ;
	int power ;
	int prime ;
	int bits ;
} facCint;

typedef struct{
	unsigned testing ;
	unsigned silent ;
	unsigned help ;
	unsigned qs_limit ;
	unsigned qs_multiplier ;
	unsigned qs_rand_seed ;
} facParams;

typedef struct {

	struct {
		void * base ;
		void * now ;
	} mem;

	facParams * params ;

	struct{
		cint cint;
		qs_sm done_up_to ;
	} trial;

	facCint * number ; // the number to factor

	cint vars[10];
	cint_sheet * calc ;

	struct {
		facCint * data ;
		unsigned index ;
	} questions;

	struct {
		facCint * data ;
		unsigned index ;
	} answers;

} facCaller;

// Quadratic sieve structures

struct qs_relation {
	qs_sm id ; // definitive relations have a non-zero id.
	cint *X;
	struct {
		qs_sm *data;
		qs_sm length ;
	} Y;
	union {
		struct {
			qs_sm *data;
			qs_sm length;
		} Z;
		struct qs_relation * next ;
	} axis ;
	// axis :
	// - "Z" field is used by definitive relations.
	// - "next" is used by data that wait to be paired, it uses a linked list instead of a "Z" field.
};

typedef struct {

	// reference to the caller
	facCaller *caller;

	// computation sheet
	cint_sheet *calc;

	// numbers that are updated
	struct {
		cint N;
		cint FACTOR;
		cint X;
		cint KEY;
		cint VALUE;
		cint CYCLE;
		cint TEMP[5];
		cint MY[5]; // available for developers, never used
	} vars;

	// polynomial vars
	struct {
		cint A;
		cint B;
		cint C;
		cint D ;
		qs_sm d_bits ;
		qs_sm offset ;
		qs_sm min ;
		qs_sm span ;
		qs_sm span_half ;
		qs_sm gray_max ;
		qs_sm curves ;
	} poly;

	// constants
	struct {
		cint kN;
		cint ONE;
		cint LARGE_PRIME;
		cint MULTIPLIER;
		cint M_HALF;
	} constants;


	// system
	struct {
		qs_sm bytes_allocated;
		void *base;
		void *now;
	} mem;

	// parameters and miscellaneous vars
	qs_md adjustor;
	qs_sm multiplier;
	qs_sm n_bits;
	qs_sm kn_bits;
	struct {
		uint8_t **positions[2];
		uint8_t *sieve;
		uint8_t *flags;
		qs_sm length;
		qs_sm length_half;
		qs_sm cache_size;
	} m;
	qs_sm iterative_list[5];
	qs_sm error_bits;
	struct{
		qs_sm value ;
	}threshold;
	unsigned rand_seed;
	qs_sm sieve_again_perms;

	// useful data sharing same length
	struct {
		struct {
			qs_sm num;
			qs_sm size;
			qs_sm sqrt_kN_mod_prime;
			qs_sm root[2];
		} *data;
		qs_sm length;
	} base;

	// useful data sharing same length
	struct {
		qs_sm *A_indexes;
		struct {
			cint B_terms;
			qs_sm *A_inv_double_value_B_terms;
			qs_sm A_over_prime_mod_prime;
			qs_sm prime_index;
			qs_md prime_squared ;
		} *data;
		struct {
			qs_sm defined;
			qs_sm subtract_one;
			qs_sm double_value;
		} values;
	} s;

	qs_sm *buffer[2]; // proportional to "length of factor base" (medium or large)

	// uniqueness trees : [ relations, cycle finder, divisors of N, ]
	struct avl_manager uniqueness[3];

	// data collection made by algorithm
	struct {
		struct qs_relation **data;
		struct {
			qs_sm now;
			qs_sm prev;
			qs_sm needs;
			qs_sm reserved;
		} length;
		qs_md large_prime;
	} relations;

	// pointers to the divisors of N are kept in a flat array
	struct {
		qs_sm processing_index;
		qs_sm total_primes;
		qs_sm length;
		cint **data;
	} divisors;

	// lanczos has its own struct
	struct {
		qs_sm safe_length ;
		struct {
			struct qs_relation *relation;
			qs_sm y_length;
		} * snapshot ;
	} lanczos;

} qs_sheet;

// Front-End factor manager
static facCint **cFactor(const cint *, facParams *);
static inline int fac_special_cases(facCaller *);
static inline int fac_trial_division(facCaller *, int);
static inline int fac_perfect_checker(facCaller *);
static inline int fac_primality_checker(facCaller *);
static inline int fac_pollard_rho_63_bits(facCaller *);
static void fac_push(facCaller *, const cint *, int, int, int);

// Math
static inline int is_prime_4669921(qs_sm);
static double log_computation(double);
static inline qs_sm multiplication_modulo(qs_md, qs_md, qs_sm);
static inline qs_sm power_modulo(qs_md, qs_md, qs_sm);
static qs_sm tonelli_shanks(qs_sm, qs_sm);
static qs_sm modular_inverse(qs_sm, qs_sm);
static inline qs_md rand_64();
static inline qs_md rand_upto(qs_md );
static inline unsigned add_rand_seed(void *);

// Cint shortcuts
static inline void simple_inline_cint(cint *, size_t, void **);
static inline void simple_dup_cint(cint *, const cint *, void **);
static inline void simple_int_to_cint(cint *, qs_md);
static inline qs_md simple_cint_to_int(const cint *);

// Avl;
static inline struct avl_node *avl_cint_inserter(void *, const void *);

// System
static inline void *mem_aligned(void *);

// Misc.
static inline int facApplyCustomParam(const char *, const char *, int, unsigned *);
static inline char *facFillParams(facParams *, int, char **);
static char *facAnswerToString(facCint **);
static inline void facDisplayProgress(const char *, double);
static inline int facSortResult(const void * , const void *);

// Quadratic sieve functions
static inline qs_sm linear_param_resolution(const double [][2], qs_sm);
static inline void qs_parametrize(qs_sheet *);
static int quadraticSieve(facCaller *);
static inline int inner_continuation_condition(qs_sheet *);
static inline int outer_continuation_condition(qs_sheet *);
static inline void preparationpart1(qs_sheet *, facCaller *);
static inline void preparationpart2(qs_sheet *);
//
static inline void preparationpart3(qs_sheet *);
static inline qs_sm preparationpart3_michel(qs_sheet *qs);
//
static inline void preparationpart4(qs_sheet *);
static inline void preparationpart5(qs_sheet *);
static inline void preparationpart6(qs_sheet *);
static inline void getStartedIteration(qs_sheet *);
static inline void iterationPart1(qs_sheet *, const cint *, cint *);
static inline void iterationPart2(qs_sheet *, const cint *, cint *);
static inline void iterationPart3(qs_sheet *, const cint *, const cint *);
static inline qs_sm iterationPart4(const qs_sheet *, qs_sm nth_curve, qs_sm **, cint *);
static inline void iterationPart5(qs_sheet *, const cint *, const cint *);
static inline void iterationPart6(qs_sheet *, const cint *, const cint *, const cint *, cint *);
static inline void iterationPart7(qs_sheet *, qs_sm, const qs_sm *);
static inline void iterationPart8(qs_sheet *, qs_sm, const qs_sm *);
static int QSRegisterFactor(qs_sheet *);
static inline void register_relations(qs_sheet *, const cint *, const cint *, const cint *);
static inline void register_relation_kind_1(qs_sheet *, const cint *, const qs_sm * const restrict [4]);
static inline void register_relation_kind_2(qs_sheet *, const cint *, const cint *, const qs_sm * const restrict [4]);
static inline void finalizationPart1(qs_sheet *, const uint64_t *);
static inline void finalizationPart2(qs_sheet *);
static inline int finalizationPart3(qs_sheet *);

// Quadratic sieve Lanczos part
static inline void lanczosMulMxNNx64(const qs_sheet *, const uint64_t *, qs_sm, uint64_t *);
static inline void lanczosMulTransMxNNx64(const qs_sheet *, const uint64_t *, uint64_t *);
static void lanczosMul64xN_Nx64(const qs_sheet *, const uint64_t *, const uint64_t *, uint64_t *, uint64_t *);
static uint64_t lanczosFindNonSingular_sub(const uint64_t *, const uint64_t *, uint64_t *, uint64_t, uint64_t *);
static void lanczosMulNx6464x64_acc(qs_sheet *, const uint64_t *, const uint64_t *, uint64_t *, uint64_t *);
static void lanczosMul64x64_64x64(const uint64_t *, const uint64_t *, uint64_t *);
static void lanczosTransposeVector(qs_sheet *, const uint64_t *, uint64_t **);
static void lanczosCombineCols(qs_sheet *, uint64_t *, uint64_t *, uint64_t *, uint64_t *);
static inline void lanczosBuildArray(qs_sheet *, uint64_t **, size_t, size_t);
static inline uint64_t *lanczosBlockWorker(qs_sheet *);
static inline uint64_t * lanczosBlock(qs_sheet *);

#endif //FAC_HEADERS
