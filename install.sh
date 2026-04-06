#!/bin/bash
# DPT2 Setup Script for macOS and Linux
# Usage: ./install.sh

set -e

echo "=========================================="
echo "DPT2 Installation Script"
echo "=========================================="
echo ""

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Check Python
echo "Checking Python installation..."
if ! command -v python3 &> /dev/null; then
    echo -e "${RED}❌ Python 3 not found${NC}"
    echo "Please install Python 3.9 or higher"
    exit 1
fi

PYTHON_VERSION=$(python3 --version | awk '{print $2}')
echo -e "${GREEN}✓ Found Python $PYTHON_VERSION${NC}"

# Create virtual environment (optional but recommended)
echo ""
echo "Creating virtual environment..."
if [ ! -d "venv" ]; then
    python3 -m venv venv
    echo -e "${GREEN}✓ Virtual environment created${NC}"
fi

# Activate virtual environment
echo "Activating virtual environment..."
source venv/bin/activate
echo -e "${GREEN}✓ Virtual environment activated${NC}"

# Install DPT2
echo ""
echo "Installing DPT2 and dependencies..."
python3 -m pip install -e .

if [ $? -eq 0 ]; then
    echo -e "${GREEN}✓ Installation successful${NC}"
else
    echo -e "${RED}❌ Installation failed${NC}"
    exit 1
fi

# Check for BLAST
echo ""
echo "Checking for BLAST+ (optional)..."
if command -v blastn &> /dev/null; then
    BLAST_VERSION=$(blastn -version | head -1)
    echo -e "${GREEN}✓ BLAST+ found: $BLAST_VERSION${NC}"
else
    echo -e "${YELLOW}⚠ BLAST+ not found (optional - needed only for --blast-screen)${NC}"
    echo "Install with: brew install blast (macOS) or sudo apt-get install ncbi-blast+ (Linux)"
fi

# Summary
echo ""
echo "=========================================="
echo -e "${GREEN}Setup Complete!${NC}"
echo "=========================================="
echo ""
echo "To activate the environment in future sessions:"
echo "  source venv/bin/activate"
echo ""
echo "To design primers:"
echo "  python scripts/design_and_score_primers.py --genes RPP30 MRGPRX1"
echo ""
echo "For more help:"
echo "  cat INSTALL.md"
echo ""
