@echo off

echo Setting up Git hooks...

REM Copy pre-commit hook to the .git/hooks directory
copy pre-commit ..\..\.git\hooks\pre-commit

echo Git hooks set up.
