#!/bin/bash

# Test script for Salmon-edgeR

echo "Testing Salmon-edgeR pipeline..."

# Test configuration validation
CONFIG="../config/pipeline_config.yaml"

if [ ! -f "$CONFIG" ]; then
    echo "✗ Config file not found"
    exit 1
fi

echo "✓ Config file exists"

# Validate YAML syntax (if python available)
if command -v python3 &> /dev/null; then
    python3 -c "import yaml; yaml.safe_load(open('$CONFIG'))" 2>/dev/null
    if [ $? -eq 0 ]; then
        echo "✓ Config file is valid YAML"
    else
        echo "✗ Invalid YAML syntax"
        exit 1
    fi
fi

echo "✓ All tests passed"
