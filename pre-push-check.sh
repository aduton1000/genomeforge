#!/bin/bash
# pre-push-check.sh — inspect what git will and won't track

echo "========================================"
echo "  FILES BEING IGNORED (by .gitignore)"
echo "========================================"
# Lists every ignored file and which rule ignores it
git status --ignored --short | grep '^!!' || echo "  (none)"

echo ""
echo "========================================"
echo "  WHICH RULE IGNORES EACH FILE"
echo "========================================"
# For every ignored path, show the exact .gitignore line responsible
git ls-files --others --ignored --exclude-standard | while read -r f; do
    printf '%-45s -> ' "$f"
    git check-ignore -v "$f"
done || echo "  (none)"

echo ""
echo "========================================"
echo "  FILES THAT WILL BE PUSHED (tracked/staged)"
echo "========================================"
git status --short

echo ""
echo "========================================"
echo "  UNTRACKED FILES (NOT ignored, will be added if you 'git add .')"
echo "========================================"
git ls-files --others --exclude-standard || echo "  (none)"
