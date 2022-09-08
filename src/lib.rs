
//! This crate is a Rust wrapper and interface to the [SPOA](https://github.com/rvaser/spoa) (simd-accelerated partial order alignment) library.
//! This library allows the efficient generation of a consensus sequence from a set of DNA or protein sequences.
//!
//! If you use this crate, please cite the original authors of SPOA:
//!
//! [Vaser, R., Sović, I., Nagarajan, N. and Šikić, M., 2017. Fast and accurate de novo genome assembly from long uncorrected reads. Genome research, 27(5), pp.737-746.](https://genome.cshlp.org/content/27/5/737)
use libc::c_char;
use std::ffi::CStr;
use std::str;

extern "C" {
    fn poa_func(
        seqs: *const *const u8,
        quals: *const *const u8,
        num_seqs: i32,
        alignment_type: i32, // 0 = local, 1 = global, 2 = gapped
        match_score: i32,
        mismatch_score: i32,
        gap_open: i32,
        gap_extend: i32,
        gap2_open: i32,
        gap2_extend: i32,
    ) -> *const c_char;
}

/// Generates a consensus sequence from a list of sequences.
/// # Arguments
///
/// * `seqs` - a vector holding the sequences (each as a null-terminated vector of u8) to form a consensus from
/// * `consensus_max_len` - The upper bound for the output consensus length. If the output consensus sequence is longer than this value, it will be truncated to this length. Setting a large value uses more memory and runtime, since a buffer of this size is allocated internally.
/// * `alignment_type` - alignment mode: 0 = local, 1 = global, 2 = gapped
/// * `match_score` - the match score for alignment
/// * `mismatch_score` - the mismatch score for alignment
/// * `gap_open` - the gap open score for alignment
/// * `gap_extend` - the gap extend score for alignment
///
/// # Returns
/// * returns the consensus of the input sequences as a vector of u8
///
/// # Examples
///
/// ```
///     use rust_spoa::poa_consensus;
///
///     fn test_dna_consensus() {
///        let mut seqs = vec![];
///        let mut quals = vec![];
///
///        // generated each string by adding small tweaks to the expected consensus "AATGCCCGTT"
///        // convert sequences to Vec<u8>
///        for seq in ["ATTGCCCGTT\0",
///            "AATGCCGTT\0",
///            "AATGCCCGAT\0",
///            "AACGCCCGTC\0",
///            "AGTGCTCGTT\0",
///            "AATGCTCGTT\0"].iter() {
///            seqs.push((*seq).bytes().map(|x|{x as u8}).collect::<Vec<u8>>());
///        }
///        for qual in ["FFFFFFFFF\0",
///            "FFFFFFFFF\0",
///            "FFFFFFFFFF\0",
///            "FFFFFFFFFF\0",
///            "FFFFFFFFFF\0",
///            "FFFFFFFFFF\0"].iter() {
///            quals.push((*qual).bytes().map(|x|{x as u8}).collect::<Vec<u8>>());
///        }
///
///        // generate consensus sequence
///        let consensus = poa_consensus(&seqs, &quals, 1, 5, -4, -3, -1, -3, -1);
///
///    }
/// ```

pub fn poa_consensus <'a>(
    seqs: &'a Vec<Vec<u8>>,
    quals: &'a Vec<Vec<u8>>,
    alignment_type: i32,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
    gap2_open: i32,
    gap2_extend: i32,
) -> &'a str {

    if seqs.len() == 0 {
        return ""
    }

    if seqs.len() != quals.len() {
        panic!("Input sequences and qualities must be of same length");
    }

    let num_seqs = seqs.len() as i32;

    let mut seq_ptrs: Vec<*const u8> = Vec::with_capacity(seqs.len());
    let mut qual_ptrs: Vec<*const u8> = Vec::with_capacity(quals.len());

    for seq in seqs {
        if !(seq[seq.len()-1] == '\0' as u8) {
            panic!("Input sequences must be null terminated");
        }
        seq_ptrs.push(seq.as_ptr());
    }
    for qual in quals {
        if !(qual[qual.len()-1] == '\0' as u8) {
            panic!("Input qualities must be null terminated");
        }
        qual_ptrs.push(qual.as_ptr());
    }

    let c_buf: *const c_char = unsafe {
        poa_func(
            seq_ptrs.as_ptr(),
            qual_ptrs.as_ptr(),
            num_seqs,
            alignment_type,
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
            gap2_open,
            gap2_extend,
        )
    };
    let c_str: &CStr = unsafe { CStr::from_ptr(c_buf) };
    let str_slice: &str = c_str.to_str().unwrap();

    str_slice
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_poa() {
        let seqs = vec!["ATTGCCCGTT",
        "AATGCCGTT",
        "AATGCCCGAT",
        "AACGCCCGTC",
        "AGTGCTCGTT",
        "AATGCTCGTT"];
        let quals = vec!["FFFFFFFFFF",
        "FFFFFFFFF",
        "FFFFFFFFFF",
        "FFFFFFFFFF",
        "FFFFFFFFFF",
        "FFFFFFFFFF",
        ];
        let mut cseqs: Vec<Vec<u8>> = Vec::with_capacity(seqs.len());
        for seq in seqs {
            let mut tmp_seq = String::from(seq);
            tmp_seq.push_str("\0");
            cseqs.push((tmp_seq.into_bytes()).to_vec());
        }

        let mut cquals: Vec<Vec<u8>> = Vec::with_capacity(quals.len());
        for qual in quals {
            let mut tmp_qual = String::from(qual);
            tmp_qual.push_str("\0");
            cquals.push((tmp_qual.into_bytes()).to_vec());
        }

        let consensus = poa_consensus(&cseqs, &cquals, 1, 5, -4, -3, -1, -3, -1);

        let expected = "AATGCCCGTT";
        assert_eq!(consensus, expected);
    }

    #[test]
    fn test_dna_consensus() {
        let mut seqs = vec![];
        let mut quals = vec![];

        // generated each string by adding small tweaks to the expected consensus "AATGCCCGTT"
        for seq in ["ATTGCCCGTT\0",
            "AATGCCGTT\0",
            "AATGCCCGAT\0",
            "AACGCCCGTC\0",
            "AGTGCTCGTT\0",
            "AATGCTCGTT\0"].iter() {
            seqs.push((*seq).bytes().map(|x|{x as u8}).collect::<Vec<u8>>());
        }
        for qual in ["FFFFFFFFFF\0",
            "FFFFFFFFF\0",
            "FFFFFFFFFF\0",
            "FFFFFFFFFF\0",
            "FFFFFFFFFF\0",
            "FFFFFFFFFF\0"].iter() {
            quals.push((*qual).bytes().map(|x|{x as u8}).collect::<Vec<u8>>());
        }

        let consensus = poa_consensus(&seqs, &quals, 1, 5, -4, -3, -1, -3, -1);

        let expected = "AATGCCCGTT";
        assert_eq!(consensus, expected);
    }


    #[test]
    fn test_protein_consensus() {
        let mut seqs = vec![];
        let mut quals = vec![];

        // expect consensus "FNLKPSWDDCQ"
        for seq in ["FNLKESWDDCQ\0".to_string(),
            "FNLKPSWDCQ\0".to_string(),
            "FNLKSPSWDDCQ\0".to_string(),
            "FNLKASWCQ\0".to_string(),
            "FLKPSWDDCQ\0".to_string(),
            "FNLKPSWDADCQ\0".to_string()].iter() {
            seqs.push(seq.chars().into_iter().map(|x|{x as u8}).collect::<Vec<u8>>());
        }

        for qual in ["FFFFFFFFFFF\0".to_string(),
            "FFFFFFFFFF\0".to_string(),
            "FFFFFFFFFFFF\0".to_string(),
            "FFFFFFFFF\0".to_string(),
            "FFFFFFFFFF\0".to_string(),
            "FFFFFFFFFFFF\0".to_string()].iter() {
            quals.push(qual.chars().into_iter().map(|x|{x as u8}).collect::<Vec<u8>>());
        }

        let consensus = poa_consensus(&seqs, &quals, 1, 5, -4, -3, -1, -3, -1);
        eprintln!("{:?}", &consensus);

        let expected = "FNLKPSWDDCQ";
        assert_eq!(consensus, expected);

    }

    #[test]
    fn test_qualitites() {
        let mut seqs = vec![];
        let mut quals = vec![];

        // expect consensus "FNLKPSWDDCQ"
        for seq in [
            "ATTGCCCATT\0".to_string(),
            "ATTGCCCGTT\0".to_string(),
            "ATTGCCCATT\0".to_string(),
            "ATTGCCCGTT\0".to_string(),
        ].iter() {
            seqs.push(seq.chars().into_iter().map(|x|{x as u8}).collect::<Vec<u8>>());
        }

        for qual in [
            "FFFFFFFIFF\0".to_string(),
            "FFFFFFF#FF\0".to_string(),
            "FFFFFFFIFF\0".to_string(),
            "FFFFFFF#FF\0".to_string(),
        ].iter() {
            quals.push(qual.chars().into_iter().map(|x|{x as u8}).collect::<Vec<u8>>());
        }

        let consensus = poa_consensus(&seqs, &quals, 1, 5, -4, -3, -1, -3, -1);
        eprintln!("{:?}", &consensus);

        let expected = "ATTGCCCATT";
        assert_eq!(consensus, expected);
    }

    #[test]
    #[should_panic]
    fn test_not_null_terminated() {
        let mut seqs = vec![];

        // generated each string by adding small tweaks to the expected consensus "AATGCCCGTT"
        for seq in ["ATTGCCCGTT",
            "AATGCCGTT",
            "AATGCCCGAT",
            "AACGCCCGTC",
            "AGTGCTCGTT",
            "AATGCTCGTT"].iter() {
            seqs.push((*seq).bytes().map(|x|{x as u8}).collect::<Vec<u8>>());
        }

        poa_consensus(&seqs, &seqs, 1, 5, -4, -3, -1, -3, -1);

    }
}
