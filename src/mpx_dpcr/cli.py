"""Command-line interface for design_and_score_primers."""
import sys
import os

# Add scripts directory to path so we can import the script
script_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'scripts')
sys.path.insert(0, script_dir)

from design_and_score_primers import main

if __name__ == "__main__":
    main()
