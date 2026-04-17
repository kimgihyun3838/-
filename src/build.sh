#!/bin/bash
echo "Building ME551 Solver..."
pip install pyinstaller
pyinstaller me551_solver.spec --clean
echo "Done! Executable is in dist/me551_solver/"
