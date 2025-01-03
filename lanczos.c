// Block Lanczos algorithm : "nullspace over a finite field", take advantage of sparsity of the matrix.

// As an undergraduate student, this project is part of my computer science + maths training.
// If you think this Lanczos Block implementation works, please star the Jason Papadopoulos's "msieve" github.

// results of operations are last AND non-const arguments.

void lanczosMulMxNNx64(const qs_sheet *qs, const uint64_t *X, const qs_sm max_size, uint64_t *Y) {
	assert(Y != X);
	memset(Y, 0, max_size * sizeof(uint64_t));
	for (qs_sm a = 0, b ; a < qs->relations.length.now; ++a) {
		struct qs_relation *const rel = qs->relations.data[a];
		for (b = 0; b < rel->Y.length; ++b)
			Y[rel->Y.data[b]] ^= X[a];
	}
}

void lanczosMulTransMxNNx64(const qs_sheet *qs, const uint64_t *X, uint64_t *Y) {
	assert(Y != X);
	for (qs_sm a = 0, b; a < qs->relations.length.now; ++a) {
		struct qs_relation *const rel = qs->relations.data[a];
		for (Y[a] = 0, b = 0; b < rel->Y.length; ++b)
			Y[a] ^= X[rel->Y.data[b]];
	}
}

void lanczosMul64xN_Nx64(const qs_sheet *qs, const uint64_t *X, const uint64_t *Y, uint64_t *Z, uint64_t *T) {
	assert(X != Z && Y != Z);
	qs_sm a, b, c, d;
	memset(Z, 0, 256 * 8 * sizeof(*Z));
	memset(T, 0, 64 * sizeof(*T));
	for(a = 0; a < qs->relations.length.now; ++a) {
		const uint64_t tmp = X[a]; // read while writing ?!
		for (b = 0, c = 0; c < 64; c += 8, b += 256)
			Z[b + (tmp >> c & 0xff)] ^= Y[a];
	}
	for (a = 0; a < 8; ++a, ++T) {
		uint64_t tmp[8] = {0};
		for (b = 0; b < 256; ++b)
			if (b >> a & 1)
				for(c = d = 0; c < 8; ++c, d += 256)
					tmp[c] ^= Z[b + d];
		for (b = 0, c = 0; b < 8; ++b, c += 8)
			T[c] = tmp[b];
	}
}

uint64_t lanczosFindNonSingular_sub(const uint64_t *t, const uint64_t * last_s, uint64_t *s, uint64_t last_dim, uint64_t *w) {
	uint64_t i, j, dim, cols[64];
	uint64_t M[64][2], mask, *row_i, *row_j, m_0, m_1;
	for (i = 0; i < 64; ++i)
		M[i][0] = t[i], M[i][1] = (uint64_t) 1 << i;
	mask = 0;
	for (i = 0; i < last_dim; ++i)
		mask |= (uint64_t) 1 << (cols[63 - i] = last_s[i]);
	for (i = j = 0; i < 64; ++i)
		if (!(mask & ((uint64_t) 1 << i)))
			cols[j++] = i;
	for (i = dim = 0; i < 64; ++i) {
		mask = (uint64_t) 1 << (cols[i]);
		row_i = M[cols[i]];
		for (j = i; j < 64; ++j) {
			row_j = M[cols[j]];
			if (row_j[0] & mask) {
				m_0 = row_j[0];
				m_1 = row_j[1];
				row_j[0] = row_i[0];
				row_j[1] = row_i[1];
				row_i[0] = m_0;
				row_i[1] = m_1;
				break;
			}
		}
		if (j < 64) {
			for (j = 0; j < 64; ++j) {
				row_j = M[cols[j]];
				if (row_i != row_j && (row_j[0] & mask))
					row_j[0] ^= row_i[0], row_j[1] ^= row_i[1];
			}
			s[dim++] = cols[i];
			continue;
		}
		for (j = i; j < 64; ++j) {
			row_j = M[cols[j]];
			if (row_j[1] & mask) {
				m_0 = row_j[0];
				m_1 = row_j[1];
				row_j[0] = row_i[0];
				row_j[1] = row_i[1];
				row_i[0] = m_0;
				row_i[1] = m_1;
				break;
			}
		}
		if (j == 64) return 0; // sub-matrix is not invertible

		for (j = 0; j < 64; ++j) {
			row_j = M[cols[j]];
			if (row_i != row_j && (row_j[1] & mask))
				row_j[0] ^= row_i[0], row_j[1] ^= row_i[1];
		}

		row_i[0] = row_i[1] = 0;
	}
	for (i = 0; i < 64; ++i)
		w[i] = M[i][1];
	mask = 0;
	for (i = 0; i < dim; ++i)
		mask |= (uint64_t) 1 << s[i];
	for (i = 0; i < last_dim; ++i)
		mask |= (uint64_t) 1 << last_s[i];
	dim *= mask == -(uint64_t) 1;
	return dim;
}

void lanczosMulNx6464x64_acc(qs_sheet *qs, const uint64_t *X, const uint64_t *Y, uint64_t *Z, uint64_t *T) {
	qs_sm a ;
	uint64_t b, c, d, e ;
	for (b = 0; b < 8; Y += 8, Z += 256, ++b)
		for (c = 0; c < 256; ++c)
			for (d = Z[c] = 0, e = c; e; e >>= 1, ++d)
				Z[c] ^= (e & 1) * Y[d];
	for (a = 0, Z -= 2048; a < qs->relations.length.now; ++a)
		for(c = d = 0; c < 64; c += 8, d += 256) {
			const uint64_t w = X[a];
			T[a] ^= Z[d + (w >> c & 0xff)];
		}
}

void lanczosMul64x64_64x64(const uint64_t *X, const uint64_t *Y, uint64_t *Z) {
	uint64_t a, b, c, d, tmp[64];
	for (a = 0; a < 64; tmp[a++] = c) {
		for (b = 0, c = 0, d = X[a]; d; d >>= 1, ++b)
			c ^= (d & 1) * Y[b];
	}
	memcpy(Z, tmp, sizeof(tmp));
}

void lanczosTransposeVector(qs_sheet *qs, const uint64_t *X, uint64_t **Y) {
	qs_sm a ; // Z will be zeroed automatically by the loop.
	uint64_t b, c, d, * Z ;
	Z = memcpy(qs->mem.now, X, qs->relations.length.now * sizeof(*X));
	for (a = 0; a < qs->relations.length.now; ++a)
		for (b = 0, c = a >> 6, d = (uint64_t) 1 << (a % 64); Z[a]; Z[a] >>= 1, ++b)
			Y[b][c] |= (Z[a] & 1) * d;
}

void lanczosCombineCols(qs_sheet *qs, uint64_t *x, uint64_t *v, uint64_t *ax, uint64_t *av) {
	qs_sm i, j, bit_pos, num_deps ;
	uint64_t k, col, col_words ;
	uint64_t mask, *mat_1[128], *mat_2[128], *tmp;
	num_deps = 64 << (v && av);
	col_words = (qs->relations.length.now + 63) / 64;
	for (i = 0; i < num_deps; ++i) {
		mat_1[i] = qs->mem.now;
		mat_2[i] = mat_1[i] + col_words;
		qs->mem.now = mat_2[i] + col_words;
	}
	lanczosTransposeVector(qs, x, mat_1);
	lanczosTransposeVector(qs, ax, mat_2);
	if (num_deps == 128) {
		lanczosTransposeVector(qs, v, mat_1 + 64);
		lanczosTransposeVector(qs, av, mat_2 + 64);
	}
	for (i = bit_pos = 0; i < num_deps && bit_pos < qs->relations.length.now; ++bit_pos) {
		mask = (uint64_t) 1 << (bit_pos % 64);
		col = bit_pos / 64;
		for (j = i; j < num_deps; ++j)
			if (mat_2[j][col] & mask) {
				tmp = mat_1[i];
				mat_1[i] = mat_1[j];
				mat_1[j] = tmp;
				tmp = mat_2[i];
				mat_2[i] = mat_2[j];
				mat_2[j] = tmp;
				break;
			}
		if (j == num_deps)
			continue;
		for (++j; j < num_deps; ++j)
			if (mat_2[j][col] & mask)
				for (k = 0; k < col_words; ++k) {
					mat_2[j][k] ^= mat_2[i][k];
					mat_1[j][k] ^= mat_1[i][k];
				}
		++i;
	}

	for (j = 0; j < qs->relations.length.now; ++j) {
		uint64_t word = 0;
		col = j / 64;
		mask = (uint64_t) 1 << (j % 64);
		for (k = i; k < 64; ++k)
			if (mat_1[k][col] & mask)
				word |= (uint64_t) 1 << k;
		x[j] = word;
	}
	char * open = (char*)mat_1[0], * close = qs->mem.now ;
	qs->mem.now = memset(open, 0, close - open);
}

void lanczosBuildArray(qs_sheet *qs, uint64_t ** target, const size_t quantity, const size_t size){
	// ensure it remains memory for linear algebra
	const char * mem_needs = (char *)qs->mem.now + quantity * size * sizeof(uint64_t);
	const char * mem_ends = (char *)qs->mem.now + qs->mem.bytes_allocated;
	assert(mem_needs < mem_ends);
	for(size_t i = 0; i < quantity; ++i)
		target[i] = qs->mem.now, qs->mem.now = mem_aligned(target[i] + size);
}

uint64_t *lanczosBlockWorker(qs_sheet *qs) {
	const qs_sm n_cols = qs->relations.length.now, v_size = n_cols > qs->base.length ? n_cols : qs->base.length;
	const uint64_t safe_size = qs->lanczos.safe_length;
	uint64_t *md[6], *xl[2], *sm[13], *tmp, *res, mask_0, mask_1, endless_guard = 1 << 11 ;
	qs_sm i, dim_0, dim_1 ;
	qs->mem.now = mem_aligned((uint64_t*) qs->mem.now + 1) ; // keep some padding.
	lanczosBuildArray(qs, md, 6, safe_size);
	lanczosBuildArray(qs, sm, 13, 64);
	lanczosBuildArray(qs, xl, 2, safe_size < 2048 ? 2048 : safe_size);
	res = (*md) - 1; // simple "trick" to return mask + null_rows
	for (i = 0; i < 64; ++i)
		sm[12][i] = i;
	dim_0 = 0;
	dim_1 = 64;
	mask_1 = (uint64_t) -1;
	for (i = 0; i < qs->relations.length.now; ++i)
		md[1][i] = (uint64_t) rand_64();
	memcpy(md[0], md[1], v_size * sizeof(uint64_t));
	lanczosMulMxNNx64(qs, md[1], v_size, xl[1]);
	lanczosMulTransMxNNx64(qs, xl[1], md[1]);
	memcpy(xl[0], md[1], v_size * sizeof(uint64_t));
	do{
		lanczosMulMxNNx64(qs, md[1], v_size, xl[1]);
		lanczosMulTransMxNNx64(qs, xl[1], md[4]);
		lanczosMul64xN_Nx64(qs, md[1], md[4], xl[1], sm[3]);
		lanczosMul64xN_Nx64(qs, md[4], md[4], xl[1], sm[5]);
		for (i = 0; i < 64 && !(sm[3][i]); ++i);
		if (i != 64) {
			dim_0 = (qs_sm) lanczosFindNonSingular_sub(sm[3], sm[12], sm[11], dim_1, sm[0]);
			if (dim_0) {
				mask_0 = 0;
				for (i = 0; i < dim_0; ++i)
					mask_0 |= (uint64_t) 1 << sm[11][i];
				for (i = 0; i < 64; ++i)
					sm[7][i] = (sm[5][i] & mask_0) ^ sm[3][i];
				lanczosMul64x64_64x64(sm[0], sm[7], sm[7]);
				for (i = 0; i < 64; ++i)
					sm[7][i] ^= (uint64_t) 1 << i;
				lanczosMul64x64_64x64(sm[1], sm[3], sm[8]);
				for (i = 0; i < 64; ++i)
					sm[8][i] &= mask_0;
				lanczosMul64x64_64x64(sm[4], sm[1], sm[9]);
				for (i = 0; i < 64; ++i)
					sm[9][i] ^= (uint64_t) 1 << i;
				lanczosMul64x64_64x64(sm[2], sm[9], sm[9]);
				for (i = 0; i < 64; ++i)
					sm[10][i] = ((sm[6][i] & mask_1) ^ sm[4][i]) & mask_0;
				lanczosMul64x64_64x64(sm[9], sm[10], sm[9]);
				for (i = 0; i < qs->relations.length.now; ++i)
					md[4][i] &= mask_0;
				lanczosMulNx6464x64_acc(qs, md[1], sm[7], xl[1], md[4]);
				lanczosMulNx6464x64_acc(qs, md[3], sm[8], xl[1], md[4]);
				lanczosMulNx6464x64_acc(qs, md[2], sm[9], xl[1], md[4]);
				lanczosMul64xN_Nx64(qs, md[1], xl[0], xl[1], sm[7]);
				lanczosMul64x64_64x64(sm[0], sm[7], sm[7]);
				lanczosMulNx6464x64_acc(qs, md[1], sm[7], xl[1], md[0]);
				tmp = md[2], md[2] = md[3], md[3] = md[1], md[1] = md[4], md[4] = tmp;
				tmp = sm[2], sm[2] = sm[1], sm[1] = sm[0], sm[0] = tmp;
				tmp = sm[4], sm[4] = sm[3], sm[3] = tmp;
				tmp = sm[6], sm[6] = sm[5], sm[5] = tmp;
				memcpy(sm[12], sm[11], 64 * sizeof(uint64_t));
				mask_1 = mask_0;
				dim_1 = dim_0;
			}
		}
	} while(dim_0 && i != 64 && --endless_guard);

	// it sometimes succeeds at 400+ iterations during software testing.
	dim_0 *= endless_guard != 0 ;

	// ===== answer finalization =====
	// result will be a simple array of the form [mask, null_rows...]
	// it's assumed that a null mask means "miss, no answer"

	*res = 0; // mask

	if(dim_0) {
		lanczosMulMxNNx64(qs, md[0], v_size, md[3]);
		lanczosMulMxNNx64(qs, md[1], v_size, md[2]);
		lanczosCombineCols(qs, md[0], md[1], md[3], md[2]);
		lanczosMulMxNNx64(qs, md[0], v_size, md[1]);
		if (*md[1] == 0) // should hold (the buffer must contain only zero)
			if (memcmp(md[1], md[1] + 1, (v_size - 1) * sizeof(uint64_t)) == 0)
				for (i = 0; i < n_cols; *res |= (*md)[i++]);
	}

	// if no mask found : clears everything, otherwise leave [mask + null rows]
	char * open = (char*) md[*res != 0], * close = qs->mem.now ;
	qs->mem.now = memset(open, 0, close - open);
	return res;
}

void lanczosReduceMatrix(qs_sheet *qs) {
	// a filtering is not always necessary to make "lanczosBlockWorker" succeed :
	// - it writes to the relations [ Y lengths, relation counter ] will change
	qs_sm a, b, c, row, col, reduced_rows = qs->base.length, *counts;
	counts = memset(qs->buffer[1], 0, qs->base.length * sizeof(*qs->buffer[1]));
	if (qs->sieve_again_perms)
		for (a = 0; a < qs->relations.length.now; ++a) {
			// "snapshot" pointers, so they can be restored if "sieve again" is fired.
			qs->lanczos.snapshot[a].relation = qs->relations.data[a] ;
			qs->lanczos.snapshot[a].y_length = qs->relations.data[a]->Y.length ;
		}
	for (a = 0; a < qs->relations.length.now; ++a)
		for (b = 0; b < qs->relations.data[a]->Y.length; ++b)
			++counts[qs->relations.data[a]->Y.data[b]];
	//
	do {
		row = reduced_rows;
		do {
			col = qs->relations.length.now;
			for (a = b = 0; a < qs->relations.length.now; ++a) {
				struct qs_relation * restrict const rel = qs->relations.data[a];
				for (c = 0; c < rel->Y.length && counts[rel->Y.data[c]] > 1; ++c);
				if (c < rel->Y.length)
					for (; rel->Y.length;)
						--counts[rel->Y.data[--rel->Y.length]];
				else
					qs->relations.data[b++] = qs->relations.data[a]; // relation is thrown.
			}
		} while (qs->relations.length.now = b, b != col);
		for (a = reduced_rows = 0; a < qs->base.length; reduced_rows += counts[a++] != 0);
		if (b = reduced_rows + 64 , qs->relations.length.now > b) { // 64 extra rows
			for (a = b; a < qs->relations.length.now; ++a)
				for (struct qs_relation * restrict const rel = qs->relations.data[a]; rel->Y.length;)
					--counts[rel->Y.data[--rel->Y.length]];
			qs->relations.length.now = b;
		}
	} while (row != reduced_rows);
}

uint64_t *lanczosBlock(qs_sheet *qs) {
	// the worker algorithm is probabilistic with high success rate
	// it is interested in the Y field of the relations (to build its matrix)
	// it receives as input the raw relations then the reduced relations
	if (qs->n_bits == 1)
		return (uint64_t *) qs->mem.now ; // nothing to solve when N = 1, return any zeroed buffer.

	uint64_t *res;
	qs_sm tries = 4, reduce_at;
	//
	if (qs->sieve_again_perms < 3) tries <<= 1 ; // 8
	if (qs->sieve_again_perms < 2) tries <<= 1 ; // 16 tries in desperation.
	reduce_at = tries >> 1 ;
	//
	if (qs->lanczos.safe_length < qs->relations.length.now) qs->lanczos.safe_length = qs->relations.length.now ;
	if (qs->lanczos.safe_length < qs->base.length) qs->lanczos.safe_length = qs->base.length ;
	qs->lanczos.safe_length += 64 - qs->lanczos.safe_length % 64 ;
	//
	do {
		if (tries == reduce_at) // 230-bit can need reduce
			lanczosReduceMatrix(qs);
		// Try to find a subset of all exponent vectors such that
		// the sum of their exponent vectors is the zero vector.
		res = lanczosBlockWorker(qs);
	} while (!*res && --tries);
	return res;
}

// Kornél Lánczos (1893 - 1974) was a Hungarian-American and later Hungarian-Irish mathematician and physicist.

// Reportedly, he sent his thesis copy to Albert Einstein, and Einstein wrote back, saying:
// "I studied your paper as far as my present overload allowed. I believe I may say this much:
// this does involve competent and original brainwork, on the basis of which a doctorate should be obtainable ...
// I gladly accept the honorable dedication."
