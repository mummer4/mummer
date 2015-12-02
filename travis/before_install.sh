#! /bin/sh

set -e

# Install yaggo if needed
if [ ! -e "$HOME/bin/yaggo" ]; then
    mkdir -p ~/bin
    curl -L -o ~/bin/yaggo https://github.com/gmarcais/yaggo/releases/download/v1.5.9/yaggo
    chmod a+rx ~/bin/yaggo
fi

# Install automake 1.14
if [ ! -e "$HOME/bin/automake" ]; then
    wget ftp://ftp.gnu.org/gnu/automake/automake-1.14.1.tar.gz
    tar zxf automake-1.14.1.tar.gz
    cd automake-1.14.1
    ./configure --prefix="$HOME"
    make
    make install
fi

# Check exec
echo "PATH: $PATH"
echo -n "yaggo: "; which yaggo
echo -n "automake: "; which automake

# Check automake version (Should not fail!)
automake --version
automake --version | \
    ruby -ne '$_ =~ /(\d+)\.(\d+)\.(\d+)/ and ($1.to_i > 1 || ($1.to_i == 1 && $2.to_i > 13)) and exit(0); puts("Automake is too old"); exit(1)'
