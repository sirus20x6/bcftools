#!/bin/bash
# Benchmark C bcftools vs Zig bcftools-zig on WGS-scale data
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="$SCRIPT_DIR/bench_data"
C_BIN="$(cd "$SCRIPT_DIR/../.." && pwd)/bcftools"
Z_BIN="$SCRIPT_DIR/../zig-out/bin/bcftools-zig"

# Resolve LD_LIBRARY_PATH for htslib
HTSLIB_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)/htslib"
if [ ! -d "$HTSLIB_DIR" ]; then
    HTSLIB_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)/htslib"
fi

# Generate data if needed
if [ ! -f "$DATA_DIR/wgs.vcf" ]; then
    echo "Generating test data..."
    bash "$SCRIPT_DIR/gen_wgs_data.sh"
fi

# Verify binaries exist
if [ ! -x "$C_BIN" ]; then
    echo "ERROR: C bcftools not found at $C_BIN"
    echo "  Build with: cd $(dirname "$C_BIN") && make"
    exit 1
fi
if [ ! -x "$Z_BIN" ]; then
    echo "ERROR: Zig bcftools-zig not found at $Z_BIN"
    echo "  Build with: cd $SCRIPT_DIR/.. && zig build -Doptimize=ReleaseFast"
    exit 1
fi

N_VARIANTS=$(grep -c -v '^#' "$DATA_DIR/wgs.vcf")
N_GENES=$(grep -c $'\tgene\t' "$DATA_DIR/wgs.gff")

echo "=========================================="
echo "  WGS-scale csq Benchmark"
echo "=========================================="
echo "Data: $N_VARIANTS variants, $N_GENES genes, 22 chromosomes"
echo "C binary:   $C_BIN"
echo "Zig binary: $Z_BIN"
echo ""

# Number of iterations for more stable timing
N_ITER=${1:-3}
echo "Iterations per tool: $N_ITER"
echo ""

# --- C benchmark ---
echo "--- C bcftools csq ---"
C_TIMES=()
for i in $(seq 1 $N_ITER); do
    T=$( { time $C_BIN csq --force -p a -f "$DATA_DIR/wgs.fa" -g "$DATA_DIR/wgs.gff" "$DATA_DIR/wgs.vcf" -o /dev/null -Ov 2>/dev/null ; } 2>&1 )
    REAL=$(echo "$T" | grep real | awk '{print $2}')
    echo "  run $i: $REAL"
    # Extract seconds from time output (handles both 0m1.234s and 1.234s formats)
    SECS=$(echo "$REAL" | sed 's/\([0-9]*\)m/\1*60+/' | sed 's/s$//' | bc -l 2>/dev/null || echo "$REAL")
    C_TIMES+=("$SECS")
done
echo ""

# --- Zig benchmark ---
echo "--- Zig bcftools-zig csq ---"
Z_TIMES=()
for i in $(seq 1 $N_ITER); do
    T=$( { time LD_LIBRARY_PATH="$HTSLIB_DIR" $Z_BIN csq --force -p a -f "$DATA_DIR/wgs.fa" -g "$DATA_DIR/wgs.gff" "$DATA_DIR/wgs.vcf" -o /dev/null 2>/dev/null ; } 2>&1 )
    REAL=$(echo "$T" | grep real | awk '{print $2}')
    echo "  run $i: $REAL"
    SECS=$(echo "$REAL" | sed 's/\([0-9]*\)m/\1*60+/' | sed 's/s$//' | bc -l 2>/dev/null || echo "$REAL")
    Z_TIMES+=("$SECS")
done
echo ""

# --- Verification ---
echo "--- Output Verification ---"
$C_BIN csq --force -p a -f "$DATA_DIR/wgs.fa" -g "$DATA_DIR/wgs.gff" "$DATA_DIR/wgs.vcf" -Ov 2>/dev/null | grep -v "^#" | cut -f1-8 | sort > /tmp/bench_c_wgs.txt
LD_LIBRARY_PATH="$HTSLIB_DIR" $Z_BIN csq --force -p a -f "$DATA_DIR/wgs.fa" -g "$DATA_DIR/wgs.gff" "$DATA_DIR/wgs.vcf" 2>/dev/null | grep -v "^#" | cut -f1-8 | sort > /tmp/bench_z_wgs.txt

C_LINES=$(wc -l < /tmp/bench_c_wgs.txt)
Z_LINES=$(wc -l < /tmp/bench_z_wgs.txt)
echo "C output:   $C_LINES lines"
echo "Zig output: $Z_LINES lines"

if diff -q /tmp/bench_c_wgs.txt /tmp/bench_z_wgs.txt > /dev/null 2>&1; then
    echo "Result: IDENTICAL"
else
    echo "Result: DIFFERS"
    echo "First 10 differences:"
    diff /tmp/bench_c_wgs.txt /tmp/bench_z_wgs.txt | head -20
fi
echo ""
echo "=========================================="
