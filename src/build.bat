@echo off
echo Building ME551 Solver...
pip install pyinstaller
python -m PyInstaller me551_solver.spec --clean
echo Done! Executable is in dist/me551_solver/
pause
