#! /bin/sh

set -e

# Install yaggo if needed
if [ ! -e "$HOME/bin/yaggo" ]; then
    mkdir -p ~/bin
    curl -L -o ~/bin/yaggo https://github.com/gmarcais/yaggo/releases/download/v1.5.9/yaggo
    chmod a+rx ~/bin/yaggo
fi
