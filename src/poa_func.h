#ifndef POA_FUNC_H
#define POA_FUNC_H

#ifdef __cplusplus
extern "C" {
#endif
unsigned poa_func(char** seqs,  // the sequences (null-terminated) to perform multiple-sequence-alignment with.
                  char** quals, // the sequences (null-terminated) to perform multiple-sequence-alignment with.
                  int num_seqs, // the number of sequences being multiply aligned
                  int l,        // alignment mode: 0 = local align, 1 = global align, 2 = semi-global
                  int m,        // score for matching bases, e.g. 5
                  int n,        // score for mismatching bases, e.g. -4
                  int g,        // gap opening penalty (must be non-positive), e.g. -3
                  int e,        // gap extension penalty (must be non-positivie), e.g. -1
                  int q,        // gap opening penalty of the second affine function (must be non-positivie), e.g. -3
                  int c,        // gap extension penalty of the second affine function (must be non-positivie), e.g. -1
                  );


#ifdef __cplusplus
}
#endif

#endif // POA_FUNC_H
