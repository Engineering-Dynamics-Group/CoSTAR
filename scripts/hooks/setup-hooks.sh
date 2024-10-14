#!/bin/sh

echo "Setting up Git hooks..."

# Copy pre-commit hook to the .git/hooks directory
cp pre-commit ../../.git/hooks/pre-commit
chmod +x .git/hooks/pre-commit

echo "Git hooks set up."
