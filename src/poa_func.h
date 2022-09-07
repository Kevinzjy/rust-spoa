#ifndef POA_FUNC_H
#define POA_FUNC_H

#ifdef __cplusplus
extern "C" {
#endif
unsigned poa_func(char** seqs,  // the sequences (null-terminated) to perform multiple-sequence-alignment with.
                  char** quals, // the sequences (null-terminated) to perform multiple-sequence-alignment with.
                  int num_seqs, // the number of sequences being multiply aligned
                  int l,        // the alignment type: 0 = local align, 1 = global align, 2 = semi-global
                  int m,        // the score to give a sequence match in alignment, e.g. 5
                  int n,        // the score to give a sequence mismatch in alignment, e.g. -4
                  int g,        // the score to give a sequence gap in alignment, e.g. -3
                  int e,        // the score to give a sequence gap extension in alignment, e.g. -1
                  int q,        // the score to give a sequence second gap in alignment, e.g. -3
                  int c,        // the score to give a sequence second gap extension in alignment, e.g. -1
                  );


#ifdef __cplusplus
}
#endif

#endif // POA_FUNC_H
