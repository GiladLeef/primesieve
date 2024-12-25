#include "avl.c"            // the trees.
#include "cint.c"           // the integers.
#include "headers.h"        // factor headers.
#include "utils.c"          // utilities and front-end.
#include "quadratic.c"      // quadratic sieve source.
#include "lanczos.c"        // quadratic sieve Lanczos.
#include "testing.c"        // quadratic sieve tests.

// Why this project uses "cint" instead of GMP?
// - Author seeks to understand the issues of 64+ bit integers.
// - The original software goal was to factor a 200-bit RSA in 30 seconds.
// - "cint" allows us to see what is sufficient to reach the goal.

static inline void facDisplayVerbose(facCint **ans);
static inline void facDisplayHelp(const char *name);

int main(int argc, char *argv[]) {
    cint N;
    facParams config = {0};
    const char *n = facFillParams(&config, argc, argv);

    if (config.testing) {
        facMiniTests(&config);
    } else if (config.help) {
        facDisplayHelp(argv[0]);
    } else if (n) {
        const int bits = 64 + 4 * (int)strlen(n);
        cintInitByString(&N, bits, n, 10); // initialize the number as a cint.
        facCint **answer = cFactor(&N, &config); // execute the routine.
        facDisplayVerbose(answer); // print answer.
        free(answer); // release answer memory.
        free(N.mem); // release number memory.
    } else {
        fputs("usage: primesieve [-h] [-s] [number]\n", stderr);
    }

    return 0;
}

static inline void facDisplayVerbose(facCint **ans) {
    for (int i = 0; i < 100; ++i) {
        putchar(' ');
    }
    putchar('\r');
    char *str = facAnswerToString(ans);
    puts(str);
    free(str);
}

static inline void facDisplayHelp(const char *name) {
    const char *str = strrchr(name, '/');
    if (!str) str = strrchr(name, '\\');
    if (!str) str = name;
    str++; // skip the '/' or '\\' character

    puts("=== [ Welcome to the factor function help ] === \n");
    printf(" - use     ./%s 123            to see the factors of 123\n", str);
    printf(" - use     ./%s -test=150      to see a one-minute 150-bit factorization test\n", str);
    printf(" - use     ./%s -limit=250     to define a limit of bit for the quadratic sieve, default to 220-bit\n", str);
    printf(" - use     ./%s -s [number]    to not see the progress of quadratic sieve\n", str);
    putchar('\n');
}
