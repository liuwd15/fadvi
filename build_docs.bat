@echo off
REM Documentation build script for FADVI (Windows)

echo Building FADVI documentation...

REM Navigate to docs directory
cd docs

REM Clean previous build
echo Cleaning previous build...
make clean

REM Build HTML documentation
echo Building HTML documentation...
make html

REM Check if build was successful
if %ERRORLEVEL% == 0 (
    echo ✅ Documentation built successfully!
    echo 📄 Open docs/_build/html/index.html in your browser
    echo 🌐 Or run: python -m http.server -d docs/_build/html 8000
) else (
    echo ❌ Documentation build failed!
    exit /b 1
)
