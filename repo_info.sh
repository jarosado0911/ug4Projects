#!/bin/bash

# Check if a path was provided
if [ -z "$1" ]; then
    echo "Usage: $0 /path/to/repo"
    exit 1
fi

REPO_PATH="$1"

# Check if the provided path exists
if [ ! -d "$REPO_PATH" ]; then
    echo "❌ Path does not exist: $REPO_PATH"
    exit 1
fi

# Move into the repo
cd "$REPO_PATH" || exit 1

# Check if it's a Git repo
if ! git rev-parse --git-dir > /dev/null 2>&1; then
    echo "❌ Not a Git repository: $REPO_PATH"
    exit 1
fi

echo "📂 Repository Information for $REPO_PATH"
echo "-----------------------------------------"

# Latest commit hash
commit_hash=$(git rev-parse HEAD)
echo "Commit Hash:        $commit_hash"

# Version number (latest tag if available)
version_tag=$(git describe --tags --abbrev=0 2>/dev/null || echo "No tag found")
echo "Version Tag:        $version_tag"

# Date of first commit (repo creation date)
creation_date=$(git log --reverse --format="%ad" --date=iso | head -n 1)
echo "Creation Date:      $creation_date"

# Date of latest commit (last modification)
last_mod_date=$(git log -1 --format="%ad" --date=iso)
echo "Last Modified:      $last_mod_date"

# Commit count
commit_count=$(git rev-list --count HEAD)
echo "Total Commits:      $commit_count"

