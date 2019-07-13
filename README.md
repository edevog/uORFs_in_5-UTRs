# uORFs_in_5-UTRs

The programs included in this repository are designed to identify potential uORFs in 5' UTRs and assign each uORF a score. The mean log score of computationally confirmed uORFs in 5' UTRs is -2.1 with a standard deviation of 3.6.

Ribo-Seq data contains ribosome counts at each position in a genome. When ribosomes bind to start codons, they move much slower than during transcription. Therefore, when ribosomes are frozen in place on the DNA strand at any point in time, there is a larger number of ribosomes located at start sites. We use this feature of Ribo-Seq data to identify potential uORFs in 5' UTRs of human DNA. Each potential uORF is assigned a score normalized by the ribosome count in its corresponding coding region. We are currently researching other relevant features of Ribo-Seq data, RNA-seq, and biological data that could be used in conjunction to inform a statistical model to further improve accuracy of the existing algorithm.
