#!/bin/bash
# Documentation build script for FADVI

echo "Building FADVI documentation..."

# Navigate to docs directory
cd docs

# Clean previous build
echo "Cleaning previous build..."
make clean

# Build HTML documentation
echo "Building HTML documentation..."
make html

# Check if build was successful
if [ $? -eq 0 ]; then
    echo "✅ Documentation built successfully!"
    echo "📄 Open docs/_build/html/index.html in your browser"
    echo "🌐 Or run: python -m http.server -d docs/_build/html 8000"
else
    echo "❌ Documentation build failed!"
    exit 1
fi
