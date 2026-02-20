#!/bin/bash
# ==========================================================================
# CLM → Julia Ralph Loop Runner
# Autonomous AI-driven porting of CLM Fortran modules to Julia
#
# Usage:
#   ./ralph_loop.sh                    # Run all remaining modules
#   ./ralph_loop.sh --tier 2           # Run only Tier 2 modules
#   ./ralph_loop.sh --resume           # Resume from last failure
#   ./ralph_loop.sh --dry-run          # Show what would be ported
#
# Requirements:
#   - claude CLI (Claude Code) installed and authenticated
#   - Julia 1.10+ with CLM.jl project set up
#   - Fortran CLM source in ../clm/ (relative to CLM.jl)
# ==========================================================================

set -euo pipefail

# --- Configuration ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"
FORTRAN_ROOT="/Users/darrieythorsson/compHydro/code/SYMFLUENCE_data/installs/clm"
LOG_FILE="$PROJECT_DIR/PORTING_LOG.md"
STATE_FILE="$PROJECT_DIR/.ralph_state"
MAX_RETRIES=3
CLAUDE_MODEL="claude-opus-4-6"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# --- Load module list ---
source "$SCRIPT_DIR/ralph_modules.sh"

# --- Parse arguments ---
TARGET_TIER=""
RESUME=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --tier)
            TARGET_TIER="$2"
            shift 2
            ;;
        --resume)
            RESUME=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# --- Initialize log ---
if [[ ! -f "$LOG_FILE" ]]; then
    cat > "$LOG_FILE" << 'LOGHEADER'
# CLM → Julia Porting Log

## Status

| Module | Tier | Status | Attempts | Notes |
|--------|------|--------|----------|-------|
LOGHEADER
fi

# --- State management ---
get_last_completed() {
    if [[ -f "$STATE_FILE" ]]; then
        cat "$STATE_FILE"
    else
        echo "0"
    fi
}

save_state() {
    echo "$1" > "$STATE_FILE"
}

# --- Prompt template ---
generate_prompt() {
    local fortran_path="$1"
    local julia_path="$2"
    local description="$3"
    local completed_modules="$4"

    cat << PROMPT
You are porting CLM (Community Land Model) from Fortran 90 to Julia.

## Current Task
Port \`${fortran_path}\` to \`${julia_path}\`
Description: ${description}

## Context
- Fortran source: ${FORTRAN_ROOT}/${fortran_path}
- Julia target: ${PROJECT_DIR}/${julia_path}
- Julia package root: ${PROJECT_DIR}
- Already-ported modules: ${completed_modules}

## Rules
1. Read the Fortran source file completely before writing any Julia code.
2. Translate ALL subroutines, functions, and type definitions. Do not skip any.
3. Preserve Fortran variable names exactly (for traceability).
4. Replace Fortran \`pointer\` arrays with Julia typed fields in structs.
5. Replace filter-based loops (\`do fc=1,num_soilc; c=filter_soilc(fc)\`) with
   mask-based iteration or simple loops that check a boolean mask.
6. For any loop over columns or patches, write it as a function suitable for
   GPU kernels later (no global state access inside the loop body).
7. Do NOT smooth discontinuities (min/max, phase change thresholds) — keep
   them identical to Fortran.
8. Add the new file to the \`include()\` list in src/CLM.jl.
9. Write a test file in test/ that validates basic functionality.
10. Run the full test suite: \`julia --project=. -e 'using Test; include("test/runtests.jl")'\`
11. Fix any errors until all tests pass.
12. Stage and commit with message: "Port ${description} from Fortran to Julia"

## Translation Patterns

### Fortran module variables → Julia struct
\`\`\`fortran
! Fortran
type, public :: temperature_type
   real(r8), pointer :: t_soisno(:,:)
end type
\`\`\`
\`\`\`julia
# Julia
Base.@kwdef mutable struct TemperatureData
    t_soisno::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
end
\`\`\`

### Fortran filter loop → Julia mask loop
\`\`\`fortran
! Fortran
do fc = 1, num_soilc
   c = filter_soilc(fc)
   t_soisno(c,1) = ...
end do
\`\`\`
\`\`\`julia
# Julia
for c in eachindex(mask_soil)
    mask_soil[c] || continue
    t_soisno[c,1] = ...
end
\`\`\`

### Fortran subroutine → Julia function
\`\`\`fortran
subroutine Foo(bounds, num_f, filter_f, temperature_inst)
   type(bounds_type), intent(in) :: bounds
\`\`\`
\`\`\`julia
function foo!(temperature::TemperatureData, mask::BitVector, bounds::UnitRange{Int})
\`\`\`

Do not modify any previously ported module unless absolutely necessary for compilation.
PROMPT
}

# --- Main loop ---
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  CLM → Julia Ralph Loop Runner${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""
echo "Project:  $PROJECT_DIR"
echo "Fortran:  $FORTRAN_ROOT"
echo "Modules:  ${#MODULES[@]}"
echo ""

completed_list=""
start_idx=0

if $RESUME; then
    start_idx=$(get_last_completed)
    echo -e "${YELLOW}Resuming from module $start_idx${NC}"
fi

total=${#MODULES[@]}
passed=0
failed=0
skipped=0

for idx in $(seq $start_idx $((total - 1))); do
    module_spec="${MODULES[$idx]}"
    IFS=':' read -r tier fortran_path julia_path description max_lines <<< "$module_spec"

    # Filter by tier if specified
    if [[ -n "$TARGET_TIER" && "$tier" != "$TARGET_TIER" ]]; then
        skipped=$((skipped + 1))
        continue
    fi

    # Check if Julia file already exists and has content
    if [[ -f "$PROJECT_DIR/$julia_path" ]]; then
        echo -e "${YELLOW}[SKIP] $description — already exists${NC}"
        completed_list="$completed_list, $julia_path"
        skipped=$((skipped + 1))
        save_state $((idx + 1))
        continue
    fi

    echo ""
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BLUE}[$((idx + 1))/$total] Tier $tier: $description${NC}"
    echo -e "${BLUE}  Fortran: $fortran_path ($max_lines lines)${NC}"
    echo -e "${BLUE}  Julia:   $julia_path${NC}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

    if $DRY_RUN; then
        echo -e "${YELLOW}  [DRY RUN] Would port this module${NC}"
        continue
    fi

    # Ensure target directory exists
    mkdir -p "$(dirname "$PROJECT_DIR/$julia_path")"

    # Generate prompt
    prompt=$(generate_prompt "$fortran_path" "$julia_path" "$description" "$completed_list")

    # Retry loop
    retry=0
    success=false

    while [[ $retry -lt $MAX_RETRIES ]]; do
        retry=$((retry + 1))
        echo -e "  ${YELLOW}Attempt $retry/$MAX_RETRIES...${NC}"

        # Write prompt to temp file (avoids ARG_MAX limits)
        prompt_file=$(mktemp)
        echo "$prompt" > "$prompt_file"

        # Run Claude with timeout (30 min max per module)
        if timeout 1800 claude -p "$(cat "$prompt_file")" \
            2>&1 | tee "$PROJECT_DIR/.ralph_last_output.log"; then

            # Verify the Julia file was created
            if [[ -f "$PROJECT_DIR/$julia_path" ]]; then
                # Run tests
                echo -e "  ${YELLOW}Running tests...${NC}"
                cd "$PROJECT_DIR"
                if julia --project=. -e 'using Test; include("test/runtests.jl")' 2>&1; then
                    echo -e "  ${GREEN}✓ PASSED${NC}"
                    success=true

                    # Commit
                    cd "$PROJECT_DIR"
                    git add -A
                    git commit -m "$(cat <<EOF
Port $description from Fortran to Julia

Tier: $tier
Source: $fortran_path ($max_lines lines)
Target: $julia_path

Co-Authored-By: Claude Opus 4.6 <noreply@anthropic.com>
EOF
)" 2>/dev/null || true

                    break
                else
                    echo -e "  ${RED}✗ Tests failed (attempt $retry)${NC}"
                fi
            else
                echo -e "  ${RED}✗ Julia file not created (attempt $retry)${NC}"
            fi
        else
            echo -e "  ${RED}✗ Claude invocation failed (attempt $retry)${NC}"
        fi

        # Clean up temp file and pause between retries
        rm -f "$prompt_file"
        sleep 2
    done
    rm -f "$prompt_file" 2>/dev/null

    # Record result
    if $success; then
        passed=$((passed + 1))
        completed_list="$completed_list, $julia_path"
        echo "| $description | $tier | ✓ PASSED | $retry | — |" >> "$LOG_FILE"
    else
        failed=$((failed + 1))
        echo "| $description | $tier | ✗ FAILED | $MAX_RETRIES | See .ralph_last_output.log |" >> "$LOG_FILE"

        # Create WIP stub so the loop can continue
        if [[ ! -f "$PROJECT_DIR/$julia_path" ]]; then
            echo "# [WIP] $description — Ralph loop failed after $MAX_RETRIES attempts" > "$PROJECT_DIR/$julia_path"
            echo "# Source: $fortran_path" >> "$PROJECT_DIR/$julia_path"
            echo "# TODO: Manual port required" >> "$PROJECT_DIR/$julia_path"
        fi

        git add -A
        git commit -m "[WIP] Failed to port $description" 2>/dev/null || true
    fi

    save_state $((idx + 1))
done

# --- Summary ---
echo ""
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}  Ralph Loop Complete${NC}"
echo -e "${BLUE}========================================${NC}"
echo -e "  ${GREEN}Passed:  $passed${NC}"
echo -e "  ${RED}Failed:  $failed${NC}"
echo -e "  ${YELLOW}Skipped: $skipped${NC}"
echo -e "  Total:   $total"
echo ""
echo "Log: $LOG_FILE"

if [[ $failed -gt 0 ]]; then
    echo ""
    echo -e "${YELLOW}Some modules failed. You can:${NC}"
    echo "  1. Fix manually and re-run: ./ralph_loop.sh --resume"
    echo "  2. Skip to next tier: ./ralph_loop.sh --tier N"
    echo "  3. Check logs: cat .ralph_last_output.log"
    exit 1
fi
