#!/bin/sh
set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color
BOLD='\033[1m'

# Helper functions
log_step() {
    printf "\n${GREEN}${BOLD}%s${NC}\n" "$1"
}

log_warning() {
    printf "${YELLOW}${BOLD}%s${NC}\n" "$1"
}

log_error() {
    printf "${RED}${BOLD}%s${NC}\n" "$1"
}

# Function to check if a command exists
check_command() {
    if ! command -v "$1" >/dev/null 2>&1; then
        log_error "$1 is required but not installed."
        case "$1" in
            deno)
                echo "Visit https://deno.com/ to install it."
                ;;
            cargo)
                echo "Visit https://rustup.rs/ to install it."
                ;;
        esac
        exit 1
    fi
}

# Create output directory with absolute path
OUTPUT_DIR="$PWD/target/reports"
log_step "Setting up output directory: $OUTPUT_DIR"
rm -rf "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

# Export as environment variable for the Deno scripts
export LOCAL_REPORT_DIR="$OUTPUT_DIR"

# Check for required tools
log_step "Checking required tools..."
check_command deno
check_command cargo

# Run formatters and checks
log_step "Running rustfmt check..."
if ! cargo fmt --check; then
    log_error "Code format issues found. Run 'cargo fmt' to fix them."
    exit 1
fi

log_step "Running Clippy and generating report..."
if ! deno run --allow-all ci/clippy_report.ts; then
    EXIT_CODE=$?
    log_error "Clippy check failed with exit code: $EXIT_CODE"
    # Continue with test execution even if Clippy fails
fi

log_step "Running tests and generating report..."
if ! deno run --allow-all ci/test_report.ts; then
    EXIT_CODE=$?
    log_error "Tests failed with exit code: $EXIT_CODE"
    # Continue to show report locations even if tests fail
fi

# Display report locations if they exist
if [ -d "$OUTPUT_DIR" ]; then
    log_step "Reports written to:"
    if [ -f "$OUTPUT_DIR/clippy-report.json" ]; then
        printf "  - Clippy report: %s\n" "$OUTPUT_DIR/clippy-report.json"
    fi
    if [ -f "$OUTPUT_DIR/test-report.json" ]; then
        printf "  - Test report: %s\n" "$OUTPUT_DIR/test-report.json"
    fi
    if [ -d "$OUTPUT_DIR/clippy-report/annotations" ] || [ -d "$OUTPUT_DIR/test-report/annotations" ]; then
        printf "  - Annotations: %s/{clippy,test}-report/annotations/*.json\n" "$OUTPUT_DIR"
    fi
fi

# Exit with error if any of the checks failed
if [ "$EXIT_CODE" != "" ]; then
    exit $EXIT_CODE
fi

log_step "All checks completed successfully!"
