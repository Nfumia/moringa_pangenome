#!/bin/bash

echo "========================================"
echo "Verifying Moringa Pangenome Environment"
echo "========================================"
echo ""

# Check conda environment
echo "Current environment: $CONDA_DEFAULT_ENV"
if [ "$CONDA_DEFAULT_ENV" != "moringa_pangenome" ] && [ "$CONDA_DEFAULT_ENV" != "moringa_pangenome_minimal" ]; then
    echo "⚠️  WARNING: Not in moringa_pangenome environment"
    echo "   Run: conda activate moringa_pangenome"
    echo ""
fi

# List of required tools
TOOLS=(hifiasm yahs bwa-mem2 samtools seqkit blastn minimap2)

echo "Checking required software:"
echo "----------------------------"

ALL_OK=true

for tool in "${TOOLS[@]}"; do
    if command -v $tool &> /dev/null; then
        VERSION=$($tool --version 2>&1 | head -1 || echo "version unknown")
        echo "✓ $tool: $VERSION"
    else
        echo "✗ $tool: NOT FOUND"
        ALL_OK=false
    fi
done

echo ""
echo "Checking optional software:"
echo "---------------------------"

OPTIONAL=(meryl busco quast genomescope2)

for tool in "${OPTIONAL[@]}"; do
    if command -v $tool &> /dev/null; then
        echo "✓ $tool: installed"
    else
        echo "○ $tool: not installed (optional)"
    fi
done

echo ""
echo "Checking Python packages:"
echo "-------------------------"

python -c "import pandas; import numpy; import matplotlib; import seaborn; print('✓ All core Python packages installed')" 2>/dev/null || echo "✗ Some Python packages missing"

echo ""
echo "========================================"

if [ "$ALL_OK" = true ]; then
    echo "✓ Installation verified successfully!"
    echo "  All required tools are available."
else
    echo "✗ Installation incomplete"
    echo "  Some required tools are missing."
    echo "  Please reinstall environment."
fi

echo "========================================"

