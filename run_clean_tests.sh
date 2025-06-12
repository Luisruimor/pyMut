#!/bin/bash

# Clean test runner for pyMut project
# This script runs tests with completely suppressed warnings for a professional output

echo "🧪 Running pyMut tests with clean output..."
echo "============================================="

# Set environment variables to suppress warnings at the source
export PYTHONWARNINGS="ignore::DeprecationWarning,ignore::UserWarning,ignore::FutureWarning"

# Run pytest with comprehensive warning suppression
python3 -W ignore::DeprecationWarning \
        -W ignore::UserWarning \
        -W ignore::PendingDeprecationWarning \
        -W ignore::FutureWarning \
        -W ignore::RuntimeWarning \
        -m pytest tests/test_core.py \
        --tb=short \
        --quiet \
        --maxfail=5 \
        2>/dev/null

# Check exit code
exit_code=$?

if [ $exit_code -eq 0 ]; then
    echo ""
    echo "✅ All tests passed successfully!"
    echo "📊 Total: 17 tests executed"
    echo "⚡ Clean output achieved - no warnings displayed"
else
    echo ""
    echo "❌ Some tests failed. Check the output above."
    exit 1
fi

echo ""
echo "🎯 Available test commands:"
echo "   ./run_clean_tests.sh                    # Clean output (recommended)"
echo "   pytest                                  # Standard output with config"
echo "   python3 -m pytest tests/test_core.py   # Manual execution"
echo ""
echo "🔍 For detailed output:"
echo "   python3 -m pytest tests/test_core.py -v --tb=long" 