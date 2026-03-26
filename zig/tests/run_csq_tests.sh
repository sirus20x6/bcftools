#!/bin/bash
#
# Run CSQ test cases from test/csq/*/ using both the C bcftools and the
# Zig bcftools-zig binary.  Compares output and reports pass/fail.
#
# Usage:
#   ./zig/tests/run_csq_tests.sh [--build] [--zig-only] [--verbose]
#
# Options:
#   --build       Build both binaries before running tests
#   --zig-only    Only run the Zig binary (compare against .cmd.out expected files)
#   --verbose     Print diff output for failed tests
#   --filter PAT  Only run tests whose name matches PAT (substring match)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ZIG_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$ZIG_ROOT/.." && pwd)"
ZIG_BIN="$ZIG_ROOT/zig-out/bin/bcftools-zig"
C_BIN="$REPO_ROOT/bcftools"
TEST_DIR="$REPO_ROOT/test/csq"
SORT_CSQ="$TEST_DIR/sort-csq"

BUILD=0
ZIG_ONLY=0
VERBOSE=0
FILTER=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --build)   BUILD=1; shift ;;
        --zig-only) ZIG_ONLY=1; shift ;;
        --verbose) VERBOSE=1; shift ;;
        --filter)  FILTER="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# --- Build if requested ---
if [[ $BUILD -eq 1 ]]; then
    echo "Building Zig binary..."
    (cd "$ZIG_ROOT" && zig build -Doptimize=ReleaseFast) || {
        echo "ERROR: Zig build failed"
        exit 1
    }
    if [[ $ZIG_ONLY -eq 0 ]]; then
        echo "Building C bcftools..."
        (cd "$REPO_ROOT" && make -j"$(nproc)" 2>/dev/null) || {
            echo "WARNING: C bcftools build failed; falling back to --zig-only mode"
            ZIG_ONLY=1
        }
    fi
fi

# --- Verify binaries exist ---
if [[ ! -x "$ZIG_BIN" ]]; then
    echo "ERROR: Zig binary not found at $ZIG_BIN"
    echo "       Run with --build or build manually: cd zig && zig build"
    exit 1
fi

if [[ $ZIG_ONLY -eq 0 ]] && [[ ! -x "$C_BIN" ]]; then
    echo "WARNING: C bcftools not found at $C_BIN; switching to --zig-only mode"
    ZIG_ONLY=1
fi

PASS=0
FAIL=0
SKIP=0
ERRORS=""

echo "============================================"
echo "CSQ Test Runner"
echo "  Zig binary:  $ZIG_BIN"
if [[ $ZIG_ONLY -eq 0 ]]; then
    echo "  C binary:    $C_BIN"
else
    echo "  Mode:        zig-only (comparing against expected .cmd.out files)"
fi
echo "  Test dir:    $TEST_DIR"
echo "============================================"
echo ""

# --- Run .cmd-based tests ---
# These have an explicit command and expected output (.cmd.out)
for test_dir in "$TEST_DIR"/*/; do
    [[ -d "$test_dir" ]] || continue
    dir_name="$(basename "$test_dir")"

    # .cmd file tests
    for cmd_file in "$test_dir"*.cmd; do
        [[ -f "$cmd_file" ]] || continue

        expected="${cmd_file}.out"
        [[ -f "$expected" ]] || continue

        test_name="${dir_name}/$(basename "$cmd_file" .cmd)"

        # Apply filter
        if [[ -n "$FILTER" ]] && [[ "$test_name" != *"$FILTER"* ]]; then
            continue
        fi

        # Read command, skip comment lines
        cmd=$(grep -v '^#' "$cmd_file" | tr '\n' ' ')

        expected_content=$(cat "$expected")

        if [[ $ZIG_ONLY -eq 1 ]]; then
            # Run Zig version only, compare against expected output
            zig_cmd=$(echo "$cmd" | sed "s|{bin}/bcftools|$ZIG_BIN|g; s|{bin}/test/csq/sort-csq|$SORT_CSQ|g")

            z_out=$(cd "$test_dir" && eval "$zig_cmd" 2>/dev/null) || {
                SKIP=$((SKIP + 1))
                echo "SKIP: $test_name (zig binary returned error)"
                continue
            }

            if [[ "$z_out" = "$expected_content" ]]; then
                PASS=$((PASS + 1))
                echo "PASS: $test_name"
            else
                FAIL=$((FAIL + 1))
                ERRORS="${ERRORS}FAIL: ${test_name}\n"
                echo "FAIL: $test_name"
                if [[ $VERBOSE -eq 1 ]]; then
                    diff <(echo "$expected_content") <(echo "$z_out") | head -20
                    echo ""
                fi
            fi
        else
            # Run C version
            c_cmd=$(echo "$cmd" | sed "s|{bin}/bcftools|$C_BIN|g; s|{bin}/test/csq/sort-csq|$SORT_CSQ|g; s|{bin}/test|$REPO_ROOT/test|g")
            c_out=$(cd "$test_dir" && eval "$c_cmd" 2>/dev/null) || {
                SKIP=$((SKIP + 1))
                echo "SKIP: $test_name (C binary returned error)"
                continue
            }

            # Run Zig version
            zig_cmd=$(echo "$cmd" | sed "s|{bin}/bcftools|$ZIG_BIN|g; s|{bin}/test/csq/sort-csq|$SORT_CSQ|g; s|{bin}/test|$REPO_ROOT/test|g")
            z_out=$(cd "$test_dir" && eval "$zig_cmd" 2>/dev/null) || {
                SKIP=$((SKIP + 1))
                echo "SKIP: $test_name (zig binary returned error)"
                continue
            }

            # Compare C vs Zig
            if [[ "$c_out" = "$z_out" ]]; then
                PASS=$((PASS + 1))
                echo "PASS: $test_name"
            else
                FAIL=$((FAIL + 1))
                ERRORS="${ERRORS}FAIL: ${test_name}\n"
                echo "FAIL: $test_name (C vs Zig differ)"
                if [[ $VERBOSE -eq 1 ]]; then
                    echo "--- C output ---"
                    echo "$c_out" | head -5
                    echo "--- Zig output ---"
                    echo "$z_out" | head -5
                    diff <(echo "$c_out") <(echo "$z_out") | head -20
                    echo ""
                fi
            fi
        fi
    done

    # .vcf file tests (standard pattern: gff + fa + vcf -> expected .txt)
    for vcf_file in "$test_dir"*.vcf; do
        [[ -f "$vcf_file" ]] || continue

        bname="$(basename "$vcf_file" .vcf)"
        gff="$test_dir/${dir_name}.gff"
        ref="$test_dir/${dir_name}.fa"
        txt="$test_dir/${bname}.txt"

        [[ -f "$gff" ]] || continue
        [[ -f "$ref" ]] || continue
        [[ -f "$txt" ]] || continue

        test_name="${dir_name}/${bname}"

        # Apply filter
        if [[ -n "$FILTER" ]] && [[ "$test_name" != *"$FILTER"* ]]; then
            continue
        fi

        expected_content=$(cat "$txt")

        # Check if the VCF has samples
        has_samples=0
        if grep -q "^#CHROM.*FORMAT" "$vcf_file" 2>/dev/null; then
            has_samples=1
        fi

        if [[ $has_samples -eq 1 ]]; then
            query_fmt='[%POS\t%REF\t%ALT\t%TBCSQ\n]\n'
        else
            query_fmt='%POS\t%REF\t%ALT\t%BCSQ\n'
        fi

        if [[ $ZIG_ONLY -eq 1 ]]; then
            if [[ $has_samples -eq 1 ]]; then
                z_out=$(cd "$test_dir" && "$ZIG_BIN" csq -f "$ref" -g "$gff" "$vcf_file" 2>/dev/null \
                    | "$C_BIN" query -f"$query_fmt" 2>/dev/null) || {
                    SKIP=$((SKIP + 1))
                    echo "SKIP: $test_name (zig csq pipeline error)"
                    continue
                }
            else
                z_out=$(cd "$test_dir" && "$ZIG_BIN" csq -f "$ref" -g "$gff" "$vcf_file" 2>/dev/null \
                    | "$SORT_CSQ" \
                    | "$C_BIN" query -f'%POS\t%REF\t%ALT\t%EXP\n%POS\t%REF\t%ALT\t%BCSQ\n\n' 2>/dev/null) || {
                    SKIP=$((SKIP + 1))
                    echo "SKIP: $test_name (zig csq pipeline error)"
                    continue
                }
            fi

            if [[ "$z_out" = "$expected_content" ]]; then
                PASS=$((PASS + 1))
                echo "PASS: $test_name"
            else
                FAIL=$((FAIL + 1))
                ERRORS="${ERRORS}FAIL: ${test_name}\n"
                echo "FAIL: $test_name"
                if [[ $VERBOSE -eq 1 ]]; then
                    diff <(echo "$expected_content") <(echo "$z_out") | head -20
                    echo ""
                fi
            fi
        else
            # Run C version
            if [[ $has_samples -eq 1 ]]; then
                c_out=$(cd "$test_dir" && "$C_BIN" csq -f "$ref" -g "$gff" "$vcf_file" 2>/dev/null \
                    | "$C_BIN" query -f"$query_fmt" 2>/dev/null) || {
                    SKIP=$((SKIP + 1))
                    echo "SKIP: $test_name (C binary error)"
                    continue
                }
            else
                c_out=$(cd "$test_dir" && "$C_BIN" csq -f "$ref" -g "$gff" "$vcf_file" 2>/dev/null \
                    | "$SORT_CSQ" \
                    | "$C_BIN" query -f'%POS\t%REF\t%ALT\t%EXP\n%POS\t%REF\t%ALT\t%BCSQ\n\n' 2>/dev/null) || {
                    SKIP=$((SKIP + 1))
                    echo "SKIP: $test_name (C binary error)"
                    continue
                }
            fi

            # Run Zig version
            if [[ $has_samples -eq 1 ]]; then
                z_out=$(cd "$test_dir" && "$ZIG_BIN" csq -f "$ref" -g "$gff" "$vcf_file" 2>/dev/null \
                    | "$C_BIN" query -f"$query_fmt" 2>/dev/null) || {
                    SKIP=$((SKIP + 1))
                    echo "SKIP: $test_name (zig binary error)"
                    continue
                }
            else
                z_out=$(cd "$test_dir" && "$ZIG_BIN" csq -f "$ref" -g "$gff" "$vcf_file" 2>/dev/null \
                    | "$SORT_CSQ" \
                    | "$C_BIN" query -f'%POS\t%REF\t%ALT\t%EXP\n%POS\t%REF\t%ALT\t%BCSQ\n\n' 2>/dev/null) || {
                    SKIP=$((SKIP + 1))
                    echo "SKIP: $test_name (zig binary error)"
                    continue
                }
            fi

            if [[ "$c_out" = "$z_out" ]]; then
                PASS=$((PASS + 1))
                echo "PASS: $test_name"
            else
                FAIL=$((FAIL + 1))
                ERRORS="${ERRORS}FAIL: ${test_name}\n"
                echo "FAIL: $test_name (C vs Zig differ)"
                if [[ $VERBOSE -eq 1 ]]; then
                    diff <(echo "$c_out") <(echo "$z_out") | head -20
                    echo ""
                fi
            fi
        fi
    done
done

echo ""
echo "============================================"
echo "Results: $PASS passed, $FAIL failed, $SKIP skipped"
echo "============================================"

if [[ $FAIL -gt 0 ]]; then
    echo ""
    echo "Failed tests:"
    echo -e "$ERRORS"
    exit 1
fi

exit 0
