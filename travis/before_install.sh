#! /bin/sh

set -e

# Install yaggo if needed
if [ ! -e "$HOME/bin/yaggo" ]; then
    mkdir -p ~/bin
    curl -L -o ~/bin/yaggo https://github.com/gmarcais/yaggo/releases/download/v1.5.9/yaggo
    chmod a+rx ~/bin/yaggo
fi

# Check automake version
automake --version
automake --version | \
    ruby -ne '$_ =~ /(\d+)\.(\d+)\.(\d+)/ and ($1.to_i > 1 || ($1.to_i == 1 && $2.to_i > 13)) and exit(0); puts("Automake is too old"); exit(1)'
