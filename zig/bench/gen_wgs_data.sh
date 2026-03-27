#!/bin/bash
# Generate synthetic WGS-scale test data for benchmarking
# Creates: bench_data/wgs.vcf, bench_data/wgs.gff, bench_data/wgs.fa

OUT_DIR="$(dirname "$0")/bench_data"
mkdir -p "$OUT_DIR"

# Parameters
N_CHROMS=22
CHROM_LEN=1000000  # 1MB per chrom (scaled down from real 250MB)
N_GENES_PER_CHROM=50
N_VARIANTS_PER_CHROM=500
N_SAMPLES=10

# Generate FASTA reference
echo "Generating reference..."
python3 -c "
import random
random.seed(42)
for c in range(1, $N_CHROMS+1):
    print(f'>chr{c}')
    seq = ''.join(random.choice('ACGT') for _ in range($CHROM_LEN))
    for i in range(0, len(seq), 80):
        print(seq[i:i+80])
" > "$OUT_DIR/wgs.fa"

# Index FASTA
samtools faidx "$OUT_DIR/wgs.fa" 2>/dev/null || {
    # If samtools not available, create a simple .fai
    echo "samtools not found, generating .fai manually..."
    python3 -c "
# Compute correct offsets for the .fai file
# Format: NAME\tLENGTH\tOFFSET\tLINEBASES\tLINEWIDTH
import os

offset = 0
with open('$OUT_DIR/wgs.fa', 'rb') as f:
    lines = []
    for c in range(1, $N_CHROMS+1):
        name = f'chr{c}'
        # Header line: >chrN\n
        header_len = len(f'>chr{c}\n'.encode())
        offset += header_len
        seq_offset = offset
        # Sequence: $CHROM_LEN bases, 80 per line
        full_lines = $CHROM_LEN // 80
        remainder = $CHROM_LEN % 80
        seq_bytes = full_lines * 81  # 80 bases + newline
        if remainder > 0:
            seq_bytes += remainder + 1
        offset += seq_bytes
        # Recompute: read file to get exact offset
    # Simpler: just scan the file
    f.seek(0)
    content = f.read()

offset = 0
i = 0
with open('$OUT_DIR/wgs.fa.fai', 'w') as out:
    data = open('$OUT_DIR/wgs.fa', 'rb').read()
    while i < len(data):
        # Parse header
        assert data[i:i+1] == b'>', f'Expected > at {i}, got {data[i:i+1]}'
        nl = data.index(b'\n', i)
        name = data[i+1:nl].decode().split()[0]
        i = nl + 1
        seq_start = i
        seq_len = 0
        first_line_len = None
        first_line_width = None
        while i < len(data) and data[i:i+1] != b'>':
            line_start = i
            nl2 = data.index(b'\n', i)
            bases = nl2 - i
            seq_len += bases
            if first_line_len is None:
                first_line_len = bases
                first_line_width = bases + 1  # including newline
            i = nl2 + 1
        out.write(f'{name}\t{seq_len}\t{seq_start}\t{first_line_len}\t{first_line_width}\n')
"
}

# Generate GFF annotation
echo "Generating GFF..."
python3 -c "
import random
random.seed(42)
print('##gff-version 3')
for c in range(1, $N_CHROMS+1):
    chr_name = f'chr{c}'
    for g in range($N_GENES_PER_CHROM):
        gene_start = random.randint(1000, $CHROM_LEN - 20000)
        gene_end = gene_start + random.randint(5000, 15000)
        strand = random.choice(['+', '-'])
        gene_id = f'gene_{c}_{g}'
        tx_id = f'tx_{c}_{g}'

        print(f'{chr_name}\t.\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\tID={gene_id};biotype=protein_coding')
        print(f'{chr_name}\t.\tmRNA\t{gene_start}\t{gene_end}\t.\t{strand}\t.\tID={tx_id};Parent={gene_id};biotype=protein_coding')

        # Generate 3-8 exons
        n_exons = random.randint(3, 8)
        exon_starts = sorted(random.sample(range(gene_start, gene_end - 200, 200), min(n_exons, (gene_end - gene_start) // 200)))

        for i, es in enumerate(exon_starts):
            ee = min(es + random.randint(100, 300), gene_end)
            phase = 0
            print(f'{chr_name}\t.\texon\t{es}\t{ee}\t.\t{strand}\t.\tParent={tx_id}')
            print(f'{chr_name}\t.\tCDS\t{es}\t{ee}\t.\t{strand}\t{phase}\tParent={tx_id}')
" > "$OUT_DIR/wgs.gff"

# Generate VCF
echo "Generating VCF..."
python3 -c "
import random
random.seed(42)
bases = 'ACGT'
# Header
print('##fileformat=VCFv4.2')
for c in range(1, $N_CHROMS+1):
    print(f'##contig=<ID=chr{c},length=$CHROM_LEN>')
print('##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">')
samples = '\t'.join(f'SAMPLE{i}' for i in range($N_SAMPLES))
print(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}')

for c in range(1, $N_CHROMS+1):
    chr_name = f'chr{c}'
    positions = sorted(random.sample(range(1, $CHROM_LEN), $N_VARIANTS_PER_CHROM))
    for pos in positions:
        ref = random.choice(bases)
        alt = random.choice([b for b in bases if b != ref])
        # Random variant type: 80% SNP, 15% del, 5% ins
        r = random.random()
        if r < 0.15:
            ref = ref + ''.join(random.choice(bases) for _ in range(random.randint(1, 5)))
        elif r < 0.20:
            alt = alt + ''.join(random.choice(bases) for _ in range(random.randint(1, 5)))

        gts = '\t'.join(random.choice(['0/0', '0/1', '1/1', '0|1', '1|0']) for _ in range($N_SAMPLES))
        print(f'{chr_name}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gts}')
" > "$OUT_DIR/wgs.vcf"

echo "Generated:"
echo "  $OUT_DIR/wgs.fa ($(wc -c < "$OUT_DIR/wgs.fa") bytes)"
echo "  $OUT_DIR/wgs.gff ($(wc -l < "$OUT_DIR/wgs.gff") lines)"
echo "  $OUT_DIR/wgs.vcf ($(grep -c -v '^#' "$OUT_DIR/wgs.vcf") variants)"
