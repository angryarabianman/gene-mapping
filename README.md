This program maps DNA sequencing reads to a reference genome (specifically E. coli) to determine where each read originates.

To run place reads into reads.fastq, because I couldn't upload on github since it is too large. 
Then run ```main.cpp```. 
For tests run ```tests.cpp```

Algorithms used:

- Seed-and-Extend Algorithm: Instead of trying to match the entire noisy read at once, the program breaks reads into small, perfect chunks called "Seeds".

- Aho-Corasick: It uses the Aho-Corasick algorithm to instantly find all locations where these seeds appear in the genome.

- Verification: Once a seed is found, the program checks the full read against that location ("Extend"). It uses Quality Scores from the data to forgive machine errors while still penalizing real mismatches.

Also, reads are split into batches, because my laptop can't handle everything at once.

<h3>REPORT:</h3>

1. Statistics:
    - Total Reads:      22720100
    - Mapped Reads:     22385723 (98.5283%)
    - Unique Mapped:    1333245 (5.95578%)
    - Multi Mapped:     21052478 (94.0442%)

2. Quality:
    - Alignment Score Threshold: 85% (using Phred penalties)

3. Genome Coverage:
    - Covered Bases:    4591389

    - Coverage %:       98.9171%
