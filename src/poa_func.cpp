#include <string.h>
#include "spoa/spoa.hpp"

extern "C" {

    // see the C header file (poa_func.h) for detailed descriptions of each argument
    const char* poa_func(char** seqs, char** quals, int num_seqs,
        int l, int m, int n, int g, int e, int q, int c) {

        if (num_seqs == 0) {
            return (unsigned) 0;
        }

        // populate the list of sequences & qualities
        std::vector<std::string> sequences;
        std::vector<std::string> qualities;
        for (int i = 0; i < num_seqs; i++){
            sequences.push_back((std::string) seqs[i]);
            qualities.push_back((std::string) quals[i]);
        }

        auto alignment_engine = spoa::AlignmentEngine::Create(static_cast<spoa::AlignmentType>(l), // Alignment mode
                                                              (int8_t) m, // match
                                                              (int8_t) n, // mismatch
                                                              (int8_t) g, // gap open
                                                              (int8_t) e, // gap extension
                                                              (int8_t) q, // second gap open
                                                              (int8_t) c // second gap extension
                                                              );

        spoa::Graph graph{};

        // add each of the real sequences (e.g. noisy sequence reads) to the graph
        // for (const auto& it: sequences) {
        //     auto alignment = alignment_engine->Align(it, graph);
        //     graph.AddAlignment(alignment, it);
        // }
        for (int i = 0; i < num_seqs; ++i) {
            const auto& it = sequences[i];
            const auto& qu = qualities[i];
            auto alignment = alignment_engine->Align(it, graph);
            graph.AddAlignment(alignment, it, qu);
        }

        // generate the consensus sequence, assign it to the allocated memory block, and return the consensus length.
        auto cns = graph.GenerateConsensus();

        // copy consensus sequence
        char *cons_str;
        cons_str = new char [cns.size() + 1];
        strcpy (cons_str, cns.c_str());
        return cons_str;
    }
}
