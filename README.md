# curation-scripts

# Requirements
POSIX Shell & Utilities
Git
Docker

# Download
git clone github.com/anastzel/curation-scripts.git $HOME/github.com/anastzel/curation-scripts
cd $HOME/github.com/anastzel/curation-scripts

# Development
cd python/python/ && make clean check all && cd ../..
cd python/ && make clean check all && cd ..

# Staging
cd python/ && make STAGING=1
