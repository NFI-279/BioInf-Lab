const CODON_TABLE = {
  'UUU': ['Phe', 'F'], 'UUC': ['Phe', 'F'], 'UUA': ['Leu', 'L'], 'UUG': ['Leu', 'L'],
  'CUU': ['Leu', 'L'], 'CUC': ['Leu', 'L'], 'CUA': ['Leu', 'L'], 'CUG': ['Leu', 'L'],
  'AUU': ['Ile', 'I'], 'AUC': ['Ile', 'I'], 'AUA': ['Ile', 'I'], 'AUG': ['Met', 'M'],
  'GUU': ['Val', 'V'], 'GUC': ['Val', 'V'], 'GUA': ['Val', 'V'], 'GUG': ['Val', 'V'],
  'UCU': ['Ser', 'S'], 'UCC': ['Ser', 'S'], 'UCA': ['Ser', 'S'], 'UCG': ['Ser', 'S'],
  'CCU': ['Pro', 'P'], 'CCC': ['Pro', 'P'], 'CCA': ['Pro', 'P'], 'CCG': ['Pro', 'P'],
  'ACU': ['Thr', 'T'], 'ACC': ['Thr', 'T'], 'ACA': ['Thr', 'T'], 'ACG': ['Thr', 'T'],
  'GCU': ['Ala', 'A'], 'GCC': ['Ala', 'A'], 'GCA': ['Ala', 'A'], 'GCG': ['Ala', 'A'],
  'UAU': ['Tyr', 'Y'], 'UAC': ['Tyr', 'Y'], 'UAA': ['STOP', '*'], 'UAG': ['STOP', '*'],
  'CAU': ['His', 'H'], 'CAC': ['His', 'H'], 'CAA': ['Gln', 'Q'], 'CAG': ['Gln', 'Q'],
  'AAU': ['Asn', 'N'], 'AAC': ['Asn', 'N'], 'AAA': ['Lys', 'K'], 'AAG': ['Lys', 'K'],
  'GAU': ['Asp', 'D'], 'GAC': ['Asp', 'D'], 'GAA': ['Glu', 'E'], 'GAG': ['Glu', 'E'],
  'UGU': ['Cys', 'C'], 'UGC': ['Cys', 'C'], 'UGA': ['STOP', '*'], 'UGG': ['Trp', 'W'],
  'CGU': ['Arg', 'R'], 'CGC': ['Arg', 'R'], 'CGA': ['Arg', 'R'], 'CGG': ['Arg', 'R'],
  'AGU': ['Ser', 'S'], 'AGC': ['Ser', 'S'], 'AGA': ['Arg', 'R'], 'AGG': ['Arg', 'R'],
  'GGU': ['Gly', 'G'], 'GGC': ['Gly', 'G'], 'GGA': ['Gly', 'G'], 'GGG': ['Gly', 'G']
};

function analyzeGene(dnaSequence) {
  // 1. Clean input & transcribe
  const dna = dnaSequence.replace(/\s+/g, '').toUpperCase();
  const rna = dna.replaceAll('T', 'U');
  
  console.log(`DNA length: ${dna.length}`);
  console.log(`\nTranscription - RNA sequence:\n${rna}`);
  
  // 2. Find first AUG (1-based index)
  const augIdx = rna.indexOf('AUG');
  if (augIdx === -1) {
    console.log("\nFirst AUG found at position: None");
    return;
  }
  
  console.log(`\nFirst AUG found at position: ${augIdx + 1}`);
  
  // 3. Read triplets from AUG until stop codon
  const codons = [];
  const aa3Letter = [];
  const aa1Letter = [];
  let stoppedAt = "End of sequence";

  for (let i = augIdx; i < rna.length - 2; i += 3) {
    const codon = rna.substring(i, i + 3);
    const info = CODON_TABLE[codon] || ['???', '?'];

    if (info[0] === 'STOP') {
      codons.push(codon); // Reading frame includes the stop codon
      stoppedAt = codon;
      break;
    }

    codons.push(codon);
    aa3Letter.push(info[0]);
    aa1Letter.push(info[1]);
  }

  // 4. Output results
  console.log(`\nReading frame:\n${codons.join('-')}`);
  console.log(`\nTranslation - 3-letter amino acid sequence:\n${aa3Letter.join('-')}`);
  console.log(`\nTranslation - 1-letter amino acid sequence:\n${aa1Letter.join('-')}`);
  console.log(`\nTranslation stopped at: ${stoppedAt}`);
}

// --- Test ---
const dnaInput = `CCGCGTAGGCTTTTGAACTGCCAAAGGTTTACGCTGATGCTCCAGGACTTACCG`;
analyzeGene(dnaInput);