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

# Install autoconf 2.69
if [ ! -e "$HOME/bin/autoconf" ]; then
    wget http://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz
    tar zxf autoconf-2.69.tar.gz
    cd autoconf-2.69
    ./configure --prefix="$HOME"
    make
    make install
fi

# Install libtool 2.4.6
if [ ! -e "$HOME/bin/libtoolize" ]; then
    wget ftp://ftp.gnu.org/gnu/libtool/libtool-2.4.6.tar.gz
    tar zxf libtool-2.4.6.tar.gz
    cd libtool-2.4.6
    ./configure --prefix="$HOME"
    make
    make install
fi

# Check exec
echo "PATH: $PATH"
echo -n "yaggo: "; which yaggo
echo -n "automake: "; which automake
echo -n "autoconf: "; which autoconf
echo -n "libtool: "; which libtool

# Check automake && autoconf version (Should not fail!)
automake --version
automake --version | \
    ruby -ne '$_ =~ /(\d+)\.(\d+)\.(\d+)/ and ($1.to_i > 1 || ($1.to_i == 1 && $2.to_i > 13)) and exit(0); puts("Automake is too old"); exit(1)'

autoconf --version
autoconf --version | \
    ruby -ne '$_ =~ /(\d+)\.(\d+)/ and ($1.to_i > 2 || ($1.to_i == 2 && $2.to_i > 68)) and exit(0); puts("Autoconf is too old"); exit(1)'

libtool --version
libtool --version | \
    ruby -ne '$_ =~ /(\d+)\.(\d+)/ and ($1.to_i > 2 || ($1.to_i == 2 && $2.to_i > 3)) and exit(0); puts("Libtool is too old"); exit(1)'
