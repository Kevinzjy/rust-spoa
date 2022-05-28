#include <string.h>
#include "spoa/spoa.hpp"

extern "C" {

    // see the C header file (poa_func.h) for detailed descriptions of each argument
    const char* poa_func(char** seqs, int num_seqs,
        int l, int m, int n, int g, int e, int q, int c) {

        if (num_seqs == 0) {
            return (unsigned) 0;
        }

        // populate the list of sequences
        std::vector<std::string> sequences;
        for (int i = 0; i < num_seqs; i++){
            sequences.push_back((std::string) seqs[i]);
        }

        auto alignment_engine = spoa::AlignmentEngine::Create(static_cast<spoa::AlignmentType>(l),
            m, n, g, e, q, c);

        spoa::Graph graph{};
        
        // add each of the real sequences (e.g. noisy sequence reads) to the graph
        for (const auto& it: sequences) {
            auto alignment = alignment_engine->Align(it, graph);
            graph.AddAlignment(alignment, it);
        }

        // generate the consensus sequence, assign it to the allocated memory block, and return the consensus length.
        auto cns = graph.GenerateConsensus();
        
        // for (int i = 0; i < l; i++){
        //     consensus[i] = cns[i];
        // }

        // sequences.clear();

        // return (unsigned) l;
        char *cons_str;
        cons_str = new char [cns.size() + 1];
        strcpy (cons_str, cns.c_str());
        return cons_str;
    }
}
